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
        dir=$(echo {output[0]} | cut -d'/' -f 1)
        dir2=$(echo {output[0]} | cut -d'/' -f 2)
        mkdir -p $dir/$dir2
        for i in {input};
        do
        j=$(echo $i | cut -d'/' -f 2 | cut -d'.' -f 1 | cut -d'0' -f 3)
        circlator fixstart $i $dir/$dir2/$dir2$j 2>> {log}; 
        done
        """

##Run progressive Mauve to compare with reference
rule pgv_mauve:
    output:
        "results/pmauve/result.png"
    input:
        og_strain = "../Akunkeei_files/fna/H3B1-04J_genomic.fna",
        new_seqs = expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    log: "logs/03-pmauve.log"
    conda: "envs/plot_region_env.yml"
    script: "03-pgvmauve.py"

##Prokka annotations
rule prokka_annot:
    output:
        expand("results/annotations/sample{i}.fna", i = ["1", "2", "3", "4"])
    input:
        protein_list = "../Akunkeei_files/faa/H3B1-04J_protein.faa",
        assemblies = expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    log: "logs/04-prokka.log"
    conda: "envs/genome_analysis_env.yml"
    threads: 8
    shell:
        """
        > {log}
        for in in {input.assemblies};
        do
        n=$(echo $in | cut -d"/" -f 3 | cut -d"." -f 1 | cut -d"i" -f 3)
        echo $n
        prokka $in --outdir results/annotations --prefix sample$n \
        --force --addgenes --genus Apilactobacillus \
        --species kunkeei --strain H3B1-04J --gram positive \
        --usegenus --protein {input.protein_list} --cpus {threads} >> {log} 2>> {log};
        done
        """

##All vs all Blast
rule all_blast:
    output:
        expand("results/blast/sample{i}.tab", i = ["1", "2", "3", "4"])
    input:
        fnas = expand("results/annotations/sample{i}.ffn", i = ["1", "2", "3", "4"]),
        og_strain = "../Akunkeei_files/cds/H3B1-04J_cds_from_genomic.fna"
    params:
        outpath = "results/blast"
    threads: 16
    log: "logs/05-blast_prokka.log"
    conda: "envs/plot_region_env.yml"
    script: "05-blast_prokka.py"

rule reformat_fna:
    output:
        expand("data/fixed_ori/sample{i}_genomic.fna", i = ["1", "2", "3", "4"])
    input:
        expand("data/fixed_ori/fixed_ori{i}.fasta", i = ["1", "2", "3", "4"])
    log: "logs/06-reformat_fna.log"
    conda: "envs/plot_region_env.yml"
    script: "06-reformat_fna.py"

rule Phase_Finder:
    output:
        per_sample=expand("results/PhaseFinder/sample{i}_genomic.tab", i = ["1", "2", "3", "4"]),
        general="results/PhaseFinder/H3B1-04J_genomic.tab"
    input:
        expand("data/fixed_ori/sample{i}_genomic.fna", i = ["1", "2", "3", "4"]),
        os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna"
    log: "logs/07-Phase_finder.log"
    conda: "envs/samtools_env.yml"
    shell:
        """
        > {log}
        bash code/07-PhaseFinder.sh {output.general} {input} >> {log} 2>> {log}
        """

##Add new rule to align reads to the genome
rule run_bwa:
    output:
        expand("bam_files/{no}.bam", no = ["01", "02", "09", "10"])
    input:
        index=os.path.expanduser("~") + "/snpseq00064/index/H3B1-04J.fna",
        reads=input_list
    log: "logs/08-reads2bam.log"
    threads: 48
    shell:
        """
        > {log}
        outdir=$(dirname -- {output[0]})
        outputs=("" {output})
        count=1
        mkdir -p $outdir
        for file in {input.reads};
        do
        outfile=$outdir/$(basename -- ${{file%.fastq.gz}}).bam
        bwa mem -x pacbio -t {threads} {input.index} $file 2>> {log} | samtools sort > $outfile 2>> {log}
        mv $outfile ${{outputs[$count]}}
        echo "Generated ${{outputs[$count]}}" >> {log}
        ((count++))
        done
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
