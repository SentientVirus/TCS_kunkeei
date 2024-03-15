import os

isol_code = {1: "02", 2: "03", 3: "08", 4: "10"}

input_list = [f"subreads/ps_405_00{key}/demultiplex.bc10{value}_BAK8A_OA--bc10{value}_BAK8A_OA.fastq.gz" for key, value in isol_code.items()]

rule all:
    input:
        expand("analysis/ps_405_00{i}/ps_405_00{i}.polished_assembly.fasta", i = ["1", "2", "3", "4"])



##First step, save all files to 001-004.fasta
rule simplify_paths:
    output:
        expand("data/00{i}.fna", i = ["1", "2", "3", "4"])
    input:
        expand("analysis/ps_405_00{i}/ps_405_00{i}.polished_assembly.fasta", i = ["1", "2", "3", "4"])
    log: "logs/00-simplify_paths.log"
    shell:
        """
        for file in {input};
        do
        output=data/$(echo $file | cut -d'/' -f 2 | cut -d'_' -f 3).fna
        cp $file $output 2>> {log}
        done
        """

##Reverse fourth file, which is the complementary strand
rule reverse_file:
    output:
        "data/rev004.fna"
    input:
        "data/004.fna"
    log: "logs/01-reverse_complement.log"
    conda: "envs/biopython_env.yml"
    script:
        "code/01-reverse_complement.py"

##Set the oriC at the right position
rule fix_ori:
    output:
        expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    input:
        expand("data/{sample}.fna", sample = ["001", "002", "003", "rev004"])
    log: "logs/02-circularize.log"
    conda: "envs/genome_analysis_env.yml"
    shell:
        """
        > {log}
        bash code/02-circularize.sh {input} {output} {log}
        """

##Run progressive Mauve to compare with reference
rule pgv_mauve:
    output:
        "results/pmauve/result.png"
    input:
        og_strain = os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna",
        new_seqs = expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    log: "logs/03x-pmauve.log"
    conda: "envs/plot_region_env.yml"
    script: "03a-pgvmauve.py"

##Prokka annotations
rule prokka_annot:
    output:
        expand("results/annotations/sample{i}.fna", i = ["1", "2", "3", "4"])
    input:
        protein_list = os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa",
        assemblies = expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    log: "logs/03-prokka.log"
    conda: "envs/genome_analysis_env.yml"
    threads: 8
    shell:
        """
        > {log}
        bash code/03-prokka.sh {threads} {input.protein_list} {input.assemblies} >> {log} 2>> {log};
        """

##All vs all Blast
rule all_blast:
    output:
        expand("results/blast/sample{i}.tab", i = ["1", "2", "3", "4"])
    input:
        fnas = expand("results/annotations/sample{i}.ffn", i = ["1", "2", "3", "4"]),
        og_strain = os.path.expanduser("~") + "/Akunkeei_files/cds/H3B1-04J_cds_from_genomic.fna"
    params:
        outpath = "results/blast"
    threads: 16
    log: "logs/04-blast_prokka.log"
    conda: "envs/plot_region_env.yml"
    script: "04-blast_prokka.py"

##Retrieve genes that are different in the re-sequenced genomes than in the reference
rule get_differences:
    output: "differences.tab"
    input: expand("results/blast/sample{i}.tab", i = ["1", "2", "3", "4"])
    threads: 1
    log: "logs/05-get_differences.log"
    conda: "envs/plot_region_env.yml"
    script: "05-parse_blast.py"


##Change FASTA line length from 60 nts to 80 nts
rule reformat_fna:
    output:
        expand("data/fixed_ori/sample{i}_genomic.fna", i = ["1", "2", "3", "4"])
    input:
        expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    log: "logs/06-reformat_fna.log"
    conda: "envs/plot_region_env.yml"
    script: "03b-reformat_fna.py"

##Check for inversions
rule Phase_Finder:
    output:
        per_sample=expand("results/PhaseFinder/sample{i}_genomic.tab", i = ["1", "2", "3", "4"]),
        general="results/PhaseFinder/H3B1-04J_genomic.tab"
    input:
        expand("data/fixed_ori/sample{i}_genomic.fna", i = ["1", "2", "3", "4"]),
        os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna"
    log: "logs/04b-Phase_finder.log"
    conda: "envs/samtools_env.yml"
    shell:
        """
        > {log}
        bash code/04b-PhaseFinder.sh {output.general} {input} >> {log} 2>> {log}
        """

##Add new rule to align reads to the genome
rule run_bwa:
    output:
        expand("bam_files/{no}.bam", no = ["01", "02", "09", "10"])
    input:
        index=os.path.expanduser("~") + "/snpseq00064/index/H3B1-04J.fna",
        reads=input_list
    log: "logs/01c-reads2bam.log"
    threads: 48
    shell:
        """
        > {log}
        bash code/01c-reads2bam.sh {input.index} {threads} {output} {input.reads} >> {log} 2>> {log}
        """

#rule synteny:
#    input:
#        ref = expand("data/fna/H3B1-04J_{record}.fna", record = ["1", "2"]),
#        tstrain = "assembly/genome.contigs.fasta" #I need to check if I have these files
#    output:
#        expand("results/synteny_comparison/H3B1-04J_{record}.crunch", record = ["1", "2"])
#    threads: 4
#        run:
#            shell("makeblastdb -in {input.ref} -dbtype nucl"),
#            shell("blastn -query {input.tstrain} -db {input.ref} \
#            -evalue 1 -task megablast -outfmt 6 > {output}')

##Use gff files from Prokka and from the Prodigal run in Circlator to perform an all vs all Blast search
