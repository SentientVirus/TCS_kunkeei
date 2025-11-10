import os

isol_code = {1: "02", 2: "03", 3: "08", 4: "10"}

input_list = [f"subreads/ps_405_00{key}/demultiplex.bc10{value}_BAK8A_OA--bc10{value}_BAK8A_OA.fastq.gz" for key, value in isol_code.items()]

#rule all:
#    input:
#        expand("analysis/ps_405_00{i}/ps_405_00{i}.polished_assembly.fasta", i = ["1", "2", "3", "4"])

rule all:
    input:
        blast = "results/blast/differences.tab",
        phase_finder = "results/PhaseFinder/H3B1-04J_genomic.tab",
        pgvmauve = "results/pmauve/Flye/result.png",
        pgmauve_NGI = "results/pmauve/NGI/result.png",
        pgamuve_combined = "results/pmauve/combined/result.png",
        bam = expand("results/bam/{no}.bam", no = ["01", "02", "09", "10"])

##First step, running Flye to get the genomes with plasmids
rule run_Flye:
    output:
        assembly = expand("assemblies/Flye/{isolate}/assembly.fasta", isolate = ["01", "02", "09", "10"]),
        stats = expand("assemblies/Flye/{isolate}/assembly_info.txt", isolate = ["01", "02", "09", "10"])
    input:
        input_list
    threads: 24
    log: "logs/01a-run_Flye.log"
    conda: "pixi_assembly/default.yml"
    shell:
        "bash code/01a-Flye_assembly.sh {input} {output.assembly} {threads} > {log} 2> {log}"

##Second step, filter out contigs with low coverage
rule filter_contigs:
    output:
        expand("assemblies/Flye/{isolate}/assembly_filtered.fasta", isolate = ["01", "02", "09", "10"])
    input:
        assembly = expand("assemblies/Flye/{isolate}/assembly.fasta", isolate = ["01", "02", "09", "10"]),
        stats = expand("assemblies/Flye/{isolate}/assembly_info.txt", isolate = ["01", "02", "09", "10"])
    threads: 2
    log: "logs/02a-filter_contigs.log"
    conda: "pixi_genome/default.yml"
    script:
        "code/02a-filter_contigs.py"

##Set the oriC at the right position
rule fix_ori:
    output:
        expand("assemblies/Flye/fixed_ori/{isolate}.fasta", isolate = ["01", "02", "09", "10"])
    input:
        assemblies = expand("assemblies/Flye/{isolate}/assembly_filtered.fasta", isolate = ["01", "02", "09", "10"]),
        start_genes = "circlator/start_genes.fna"
    log: "logs/03a-circularize.log"
    conda: "pixi_genome/genome_analysis.yml"
    shell:
        """
        mkdir -p $(basename -- {output[0]})
        > {log}
        bash code/03-circularize.sh {input.assemblies} {output} {input.start_genes} {log}
        """

##Run progressive Mauve to compare with reference
rule pgv_mauve:
    output:
        "results/pmauve/Flye/result.png"
    input:
        og_strain = os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna",
        new_seqs = expand("assemblies/Flye/fixed_ori/{isolate}.fasta", isolate = ["01", "02", "09", "10"])
    log: "logs/04a-pmauve.log"
    conda: "pixi_genome/genome_analysis.yml"
    script: "code/04-pgvmauve.py"


##Assembly from the sequencing facility
##First step, save all files to 001-004.fasta
rule simplify_paths:
    output:
        expand("assemblies/NGI/00{i}.fasta", i = ["1", "2", "3", "4"])
    input:
        expand("analysis/ps_405_00{i}/ps_405_00{i}.polished_assembly.fasta", i = ["1", "2", "3", "4"])
    log: "logs/01b-simplify_paths.log"
    shell:
        """
        > {log}
        outdir=$(dirname -- {output[0]})
        echo "Create directory "$outdir >> {log}
        for file in {input};
        do
        echo "Path to input "$file >> {log}
        output=$outdir/$(echo $file | cut -d'/' -f 2 | cut -d'_' -f 3).fasta
        cp $file $output
        echo "Created output "$output >> {log}
        done
        """

##Reverse fourth file, which is the complementary strand
rule reverse_file:
    output:
        "assemblies/NGI/rev004.fasta"
    input:
        "assemblies/NGI/004.fasta"
    log: "logs/02b-reverse_complement.log"
    conda: "pixi_genome/default.yml"
    script:
        "code/02b-reverse_complement.py"


##Set the oriC at the right position
rule fix_ori_seq:
    output:
        expand("assemblies/NGI/fixed_ori/{i}.fasta", i = ["01", "02", "09", "10"])
    input:
        assemblies = expand("assemblies/NGI/{i}.fasta", i = ["001", "002", "003", "rev004"]),
        start_genes = "circlator/start_genes.fna"
    log: "logs/03b-circularize.log"
    conda: "pixi_genome/genome_analysis.yml"
    shell:
        """
        mkdir -p $(basename -- {output[0]})
        > {log}
        bash code/03-circularize.sh {input.assemblies} {output} {input.start_genes} {log}
        """

##Run progressive Mauve to compare with reference
rule pgv_mauve_seq:
    output:
        "results/pmauve/NGI/result.png"
    input:
        og_strain = os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna",
        new_seqs = expand("assemblies/NGI/fixed_ori/{i}.fasta", i = ["01", "02", "09", "10"])
    log: "logs/04b-pmauve.log"
    conda: "pixi_genome/genome_analysis.yml"
    script: "code/04-pgvmauve.py"


##Based on the progressiveMauve results, align the reads back to the original assembly to check if there are any true differences
rule pacbio2ref:
    output:
        index="index/H3B1-04J.fna",
        bamfiles=expand("results/bam/{no}.bam", no = ["01", "02", "09", "10"]),
        bam_index=expand("results/bam/{no}.bam.bai", no = ["01", "02", "09", "10"])
    input:
        ref=os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna",
        reads=input_list
    log: "logs/05-reads2bam.log"
    conda: "pixi_genome/genome_analysis.yml"
    threads: 48
    shell:
        """
        > {log}
        bash code/05-reads2bam.sh {input.ref} {output.index} {threads} {output.bamfiles} {input.reads} >> {log} 2>> {log}
        """

##Create combined assemblies
rule combine_assemblies:
    output:
        assemblies = expand("assemblies/combined/{isolate}.fasta", isolate = ["01", "02", "09", "10"]),
        pgvmauve = "results/pmauve/combined/result.png"
    input:
        NGI = expand("assemblies/NGI/fixed_ori/{i}.fasta", i = ["01", "02", "09", "10"]),
        Flye = expand("assemblies/Flye/fixed_ori/{i}.fasta", i = ["01", "02", "09", "10"]),
        og_strain = os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna"
    log: "logs/06-combine_assemblies.log"
    conda: "pixi_genome/genome_analysis.yml"
    threads: 1
    script: "code/06-combine_assemblies.py"

##Prokka annotations
rule prokka_annot:
    output:
        expand("results/annotations/prokka/{i}.{ext}", i = ["01", "02", "09", "10"], ext = ["err", "gbk", "gff", "faa", "fna", "ffn", "fsa", "log", "sqn", "tbl", "tsv", "txt"])
    input:
        protein_list = os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa",
        assemblies = expand("assemblies/combined/{i}.fasta", i = ["01", "02", "09", "10"])
    log: "logs/07-prokka.log"
    conda: "pixi_genome/genome_analysis.yml"
    threads: 8
    shell:
        """
        > {log}
        bash code/07-prokka.sh {threads} {input.protein_list} {input.assemblies} >> {log} 2>> {log};
        """

##All vs all Blast
rule all_blast:
    output:
        expand("results/blast/{i}.tab", i = ["01", "02", "09", "10"])
    input:
        fnas = expand("results/annotations/prokka/{i}.ffn", i = ["01", "02", "09", "10"]),
        og_strain = os.path.expanduser("~") + "/Akunkeei_files/cds/H3B1-04J_cds_from_genomic.fna"
    params:
        outpath = "results/blast"
    threads: 16
    log: "logs/08-blast_prokka.log"
    conda: "pixi_genome/genome_analysis.yml"
    script: "code/08-blast_prokka.py"

##Retrieve genes that are different in the re-sequenced genomes than in the reference
rule get_differences:
    output: "results/blast/differences.tab"
    input: expand("results/blast/{i}.tab", i = ["01", "02", "09", "10"])
    threads: 1
    log: "logs/09-get_differences.log"
    conda: "pixi_genome/genome_analysis.yml"
    script: "code/09-parse_blast.py"


##Change FASTA line length from 60 nts to 80 nts
rule reformat_fna:
    output:
        expand("assemblies/combined/{i}_80nts.fasta", i = ["01", "02", "09", "10"])
    input:
        expand("assemblies/combined/{i}.fasta", i = ["01", "02", "09", "10"])
    log: "logs/07c-reformat_fna.log"
    conda: "pixi_genome/default.yml"
    script: "code/07c-reformat_fna.py"

##Check for inversions
rule Phase_Finder:
    output:
        per_sample=expand("results/PhaseFinder/sample{i}_genomic.tab", i = ["01", "02", "09", "10"]),
        general="results/PhaseFinder/H3B1-04J_genomic.tab"
    input:
        expand("assemblies/combined/{i}_80nts.fasta", i = ["01", "02", "09", "10"]),
        os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna"
    log: "logs/08c-Phase_finder.log"
#    conda: "envs/samtools_env.yml"
    shell:
        """
        > {log}
        bash code/08c-PhaseFinder.sh {output.general} {input} >> {log} 2>> {log}
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
