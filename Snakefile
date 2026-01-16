configfile: "config.yml"
import os

isol_code = {"01-F-1": "S5", "01-F-2": "S13", "01-F-3": "S21", "01-F-4": "S29", "01-F-5": "S37", "01-S-1": "S1", "01-S-2": "S9", 
    "01-S-3": "S17", "01-S-4": "S25", "01-S-5": "S33", "02-F-1": "S6", "02-F-2": "S14", "02-F-3": "S22", "02-F-4": "S30", "02-F-5": "S38",
    "02-S-1": "S2", "02-S-2": "S10", "02-S-3": "S18", "02-S-4": "S26", "02-S-5": "S34", "09-F-1": "S7", "09-F-2": "S15", "09-F-3": "S23", 
    "09-F-4": "S31", "09-F-5": "S39", "09-S-1": "S3", "09-S-2": "S11", "09-S-3": "S19", "09-S-4": "S27", "09-S-5": "S35", "10-F-1": "S8",
    "10-F-2": "S16", "10-F-3": "S24", "10-F-4": "S32", "10-F-5": "S40", "10-S-1": "S4", "10-S-2": "S12", "10-S-3": "S20", "10-S-4": "S28",
    "10-S-5": "S36"}

input_base = []
all_input = [f"VF-3336-H3B1-04J-{key}_{value}_L001" for key, value in isol_code.items()]
summary_input = [f"I{key.replace('-', '_')[:-2]}_{value}" for key, value in isol_code.items()]


##Functions to retrieve inputs
def add_path_extension(lst, path = "", extension = "", extra = ""):
    new_list = []
    for element in lst:
        new_element = path + "/" + element + extra + "." + extension
        new_list.append(new_element)
    return new_list

def input_def(lst, path = "", extension = ""):
    new_list = []
    for element in lst:
        new_element = path + "/" + "Sample_" + element.split("_")[0] + "/" + element + "_R1_001" + "." + extension
        new_list.append(new_element)
        new_list.append(new_element.replace("R1", "R2"))
    return new_list

def input_strand(lst, strand = 1, path = "", extension = ""):
    new_list = []
    for element in lst:
        if strand == 1:
            new_element = path + "/" + "Sample_" + element.split("_")[0] + "/" + element + "_R1_001" + "." + extension
        elif strand == -1:
            new_element = path + "/" + "Sample_" + element.split("_")[0] + "/" + element + "_R2_001" + "." + extension
        new_list.append(new_element)
    return new_list


##Rule to define all desired outputs
rule all:
    input:
        multiqc = "results/multiqc/multiqc_report.html",
        improved_annot = expand("results/DE/{comparison}{ext}_improved_annot.tsv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["", "_lfc1"]),
        coverage = add_path_extension(all_input, "results/coverage", "perbase.cov"),
        picard = add_path_extension(all_input, "results/picard", "pdf", "_insert_size_histogram"),
        saturation = expand("plots/saturation{extra}_k0.ps", extra = ["_collapsed", ""]),
        volcano = expand("plots/{comparison}_volcano.{ext}", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["png", "pdf"])

##Rule to index the reference genome of strain H3B1-04J
rule index_genome:
    output:
        fna = "index/H3B1-04J.fna",
        index = expand("index/H3B1-04J.fna.{ext}", ext = ["amb", "ann", "bwt", "pac", "sa"])
    input:
        fna = os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna"
    conda: "pixi_transcript/rnaseq.yml"
    log: "logs/01-index_genome.log"
    shell:
        "bash code/01-index_genome.sh {input} {output.fna} 2> {log}"

##Rule to trim the Illumina MiSeq RNA reads from the resequenced isolates of H3B1-04J
rule trim_reads:
    output:
        R1 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair1"),
        R2 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair2")
    input:
        R1 = input_strand(all_input, strand = 1, path = "files/VF-3336/221006_M06455_0144_000000000-KMH8C", extension = "fastq.gz"),
        R2 = input_strand(all_input, strand = -1, path = "files/VF-3336/221006_M06455_0144_000000000-KMH8C", extension = "fastq.gz"),
        adapter = "adapters/adapters.fasta"
    params: outdir = "trimmed_reads"
    conda: "pixi_transcript/rnaseq.yml"
    log: add_path_extension(all_input, path = "logs/02-read_trimming", extension = "log", extra = "-trimmed")
    shell:
        """
        bash code/02-read_trimming.sh {params.outdir} {input.adapter} {input.R1} &&
        mv trimmed_reads/*.log logs/02-read_trimming
        """

##Rule to run FastQC and MultiQC on trimmed reads
rule RNA_read_quality_control:
    output:
        R1 = add_path_extension(all_input, path = "results/fastqc", extension = "zip", extra = "-trimmed-pair1_fastqc"),
        R2 = add_path_extension(all_input, path = "results/fastqc", extension = "zip", extra = "-trimmed-pair2_fastqc"),
        multiqc = "results/multiqc/multiqc_report.html"
    input: 
        R1 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair1"),
        R2 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair2")
    threads: 2
    conda: "pixi_transcript/rnaseq.yml"
    log: "logs/02.1-read_QC.log"
    shell:
        """
        bash code/02.1-read_QC.sh {threads} {output.R1[0]} {output.multiqc} {input} >> {log} 2>> {log}
        """

##Rule to align the RNA reads to the reference genome of H3B1-04J
rule align2fna:
    output:
        outbam = add_path_extension(all_input, "results/bam", "bam")
    input:
        R1 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair1"),
        R2 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair2"),
        genome = "index/H3B1-04J.fna"
    params: 
        indir = "trimmed_reads",
        outdir = "results/bam"
    log: "logs/03-read_alignment.log"
    conda: "pixi_transcript/rnaseq.yml"
    shell:
        "bash code/03-read_alignment.sh {params.outdir} {input.genome} {log} {input.R1}"

##Rule to create a gff from the new GenBank
rule gbk2gff:
    output:
        "results/gff/reference_loctag.gff"
    input:
        "../ugc00027/results/annotations/emapper2gbk/reference_loctag.gbk"
    log: "logs/04-gbk2gff.log"
    conda: "pixi_transcript/gbk2gff.yml"
    shell:
        "genbank_to -g {input} --gff3 {output} 2> {log}"
#    script: "code/04-gbk2gff.yml"


##Rule to calculate gene counts using the read alignment and the annotation of the reference strain
rule count_genes:
    output:
        counts_reverse = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        summary_reverse = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts.summary"),
        counts_forward = add_path_extension(all_input, "featureCounts_forward/nofilter", "featureCounts"),
        summary_forward = add_path_extension(all_input, "featureCounts_forward/nofilter", "featureCounts.summary")
    input:
        bam = add_path_extension(all_input, "results/bam", "bam"),
        gff = "results/gff/reference_loctag.gff" #"/Akunkeei_files/gff/H3B1-04J_genomic.gff"
    params: out1 = "featureCounts_reverse/nofilter", out2 = "featureCounts_forward/nofilter"
    log: "logs/05-read_counts.log"
    conda: "pixi_transcript/rnaseq.yml"
    shell:
        "bash code/05-read_counts.sh {params.out1} {params.out2} {input.gff} {log} {input.bam}"

##Rule to calculate RNA read coverage along the genome
rule calculate_coverage:
    output:
        add_path_extension(all_input, "results/coverage", "perbase.cov")
    input:
        add_path_extension(all_input, "results/bam", "bam")
    params: "results/coverage"
    conda: "pixi_transcript/rnaseq.yml"
    log: "logs/06-coverage.log"
    shell:
        "bash code/06-coverage.sh {params} {input} 1>&2 2> {log}"

##Rule to calculate several metrics related to the alignment
rule picard_tools:
    output:
        pdf = add_path_extension(all_input, "results/picard", "pdf", "_insert_size_histogram"),
        txt = add_path_extension(all_input, "results/picard", "txt", "_insert_size_metrics")
    input:
        add_path_extension(all_input, "results/bam", "bam")
    params: "results/picard"
    conda: "pixi_transcript/rnaseq.yml"
    log: "logs/06-picard.log"
    shell:
        "bash code/06-picard.sh {params} {input} 1>&2 2> {log}"

##Rule to calculate TPM per sample and per isolate
#rule get_TPM:
#    output:
#        out = add_path_extension(all_input, "results/TPM", "out", "_genes"),
#        ent = add_path_extension(all_input, "results/TPM", "ent", "_genes"),
#        uni = add_path_extension(all_input, "results/TPM", "uni", "_genes") 
#    input:
#        bam = add_path_extension(all_input, "results/bam", "bam"),
#        gtf = "H3B1-04J_genomic.gtf" #Modified gtf were CDS features have been renamed to exon
#    conda: "envs/circular.yml"
#    log: "logs/07-calculate_TPM.log"
#    shell:
#        "bash code/07-calculate_TPM.sh {input.gtf} {input.bam} 2> {log}"

rule get_TPM_formula:
    output:
        per_sample = "results/TPM/TPM_per_sample.tsv",
        mean = "results/TPM/mean_TPM.tsv"
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts")
    params: os.getcwd()
    log: "logs/07-calculate_TPM.log"
    conda: "pixi_transcript/default.yml"
    script:
        "code/07-calculate_TPM.py"

##Rule to prefilter out genes with low counts before the differential expression analysis
rule filter_counts:
    output:
        counts = add_path_extension(all_input, "featureCounts_reverse/filtered", "featureCounts"),
        summary = add_path_extension(summary_input, "results/summary", "tsv", "_count_distribution")
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        gbk = "../ugc00027/results/annotations/emapper2gbk/reference_loctag.gbk"  #"/Akunkeei_files/gbff/H3B1-04J_genomic.gbff" 
    conda: "pixi_transcript/default.yml"
    params: workdir = os.getcwd()
    log: "logs/08-prefilter_counts.log"
    script:
        "code/08-prefilter_counts.py"

##Rule to implement the filtering and change the count file format to a more readable format
rule parse_counts:
    output:
        avg_nofilter = "featureCounts_reverse/countfiles/nofilter/avg_gene_counts.tsv",
        avg_filtered = "featureCounts_reverse/countfiles/filtered/avg_gene_counts.tsv",
        counts = expand("featureCounts_reverse/countfiles/{dir}/H3B1-04J_{isol}{cond}_counts.tsv", dir = ["filtered", "nofilter"], isol = ["01", "02", "09", "10"], cond = ["F", "S"]),
        meta = "featureCounts_reverse/countfiles/H3B1-04J_metadata.tsv",
        pos = "meta/gene_positions.tsv"
    input:
        TPM = "results/TPM/mean_TPM.tsv",
        meta = "meta/metadata.VF-3336.csv",
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        counts_filtered = add_path_extension(all_input, "featureCounts_reverse/filtered", "featureCounts")
    params: countdir = "featureCounts_reverse/nofilter", filtered_countdir = "featureCounts_reverse/filtered"
    conda: "pixi_transcript/default.yml"
    log: "logs/09-filter_counts.log"
    script:
        "code/09-filter_counts.py"

##Rule to run a saturation analysis before the differential expression analysis
rule saturation:
    output:
        plots = expand("plots/saturation{extra}_k0.png", extra = ["_collapsed", ""]),
        ps = expand("plots/saturation{extra}_k0.ps", extra = ["_collapsed", ""])
    input:
        counts = expand("featureCounts_reverse/countfiles/filtered/H3B1-04J_{isol}{cond}_counts.tsv", isol = ["01", "02", "09", "10"], cond = ["F", "S"]),
        meta = "featureCounts_reverse/countfiles/H3B1-04J_metadata.tsv",
        pos = "meta/gene_positions.tsv"
    conda: "pixi_transcript/renv.yml"
    log: "logs/10-saturation_analysis.log"
    script:
        "code/10-saturation_analysis.R"

##Rule to run a differential expression analysis, comparing several conditions
rule differential_expression:   
    output:
        dif_expr = expand("results/DE/{comparison}{ext}.csv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["", "_lfc1"]),
        DE_annot = expand("results/DE/{comparison}{ext}_annotated.tsv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["", "_lfc1"]),
        pca = expand("plots/pcaplot.{ext}", ext = ["png", "pdf"]), dist = expand("plots/dist_plot.{ext}", ext = ["png", "ps"]),
        heatmap = expand("plots/{comparison}_heatmap.{ext}", comparison = ["global", "Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["png", "ps"]), 
        plots = expand("plots/{comparison}.{ext}", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["png", "ps"])
    input:
        counts = expand("featureCounts_reverse/countfiles/filtered/H3B1-04J_{isol}{cond}_counts.tsv", isol = ["01", "02", "09", "10"], cond = ["F", "S"]),
        meta = "featureCounts_reverse/countfiles/H3B1-04J_metadata.tsv"
    conda: "pixi_transcript/renv.yml"
    log: "logs/11-dif_expression.log"
    script:
        "code/11-dif_expression.R"

##Rule to add improved annotations to the results
rule add_improved_annotations:
    output:
        improved_annot = expand("results/DE/{comparison}{ext}_improved_annot.tsv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["", "_lfc1"])    
    input:
        ref_gbk = os.path.expanduser("~") + "/mucoid_project/ugc00027/results/annotations/emapper2gbk/reference_loctag.gbk",
        DE_annot = expand("results/DE/{comparison}{ext}_annotated.tsv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["", "_lfc1"])
    conda: "pixi_transcript/default.yml"
    log: "logs/12-modify_annot.log"
    script:
        "code/12-modify_annot.py"


##Rule to generate volcano plots
rule volcano_plots:
    output:
        volcano = expand("plots/{comparison}_volcano.{ext}", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["png", "pdf"])
    input:
        annotated_expr = expand("results/DE/{comparison}_improved_annot.tsv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"])
    conda: "pixi_transcript/renv.yml"
    log: "logs/13-volcano_plots.log"
    script: "code/13-volcano_plots.R"

