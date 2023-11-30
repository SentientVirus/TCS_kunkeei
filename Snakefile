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


rule index_genome:
    output:
        fna = "index/H3B1-04J.fna",
        index = expand("index/H3B1-04J.fna.{ext}", ext = ["amb", "ann", "bwt", "pac", "sa"]),
        ht2 = "index/H3B1-04J.{no}.ht2", no = list(range(1,9)))
    input:
        fna = os.path.expanduser("~") + "/Akunkeei_files/fna/H3B1-04J_genomic.fna"
    conda: "alignment.yml"
    log: "logs/01-index_genome.log"
    shell:
        "bash code/01-index_genome.sh {input} {output.fna} 2> {log}"

rule trim_reads:
    output:
        R1 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair1"),
        R2 = add_path_extension(all_input, path = "trimmed_reads", extension = "fastq", extra = "-trimmed-pair2")
    input:
        R1 = input_strand(all_input, strand = 1, path = "files/VF-3336/221006_M06455_0144_000000000-KMH8C", extension = "fastq.gz"),
        R2 = input_strand(all_input, strand = -1, path = "files/VF-3336/221006_M06455_0144_000000000-KMH8C", extension = "fastq.gz"),
        adapter = "TruSeq3-PE-2.fa"
    params: outdir = "trimmed_reads"
    conda: "alignment.yml"
    log: add_path_extension(all_input, path = "logs/02-read_trimming", extension = "log", extra = "-trimmed")
    shell:
        """
        bash code/02-read_trimming.sh {params.outdir} {input.adapter} {input.R1}
        mv trimmed_reads/*.log logs/02-read_trimming
        """
        
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
    conda: "alignment.yml"
    shell:
        "bash code/03-read_alignment.sh {params.outdir} {input.genome} {log} {input.R1}"

rule count_genes:
    output:
        counts_reverse = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        summary_reverse = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts.summary"),
        counts_forward = add_path_extension(all_input, "featureCounts_forward/nofilter", "featureCounts"),
        summary_forward = add_path_extension(all_input, "featureCounts_forward/nofilter", "featureCounts.summary")
    input:
        bam = add_path_extension(all_input, "results/bam", "bam"),
        gff = os.path.expanduser("~") + "/Akunkeei_files/gff/H3B1-04J_genomic.gff"
    params: out1 = "featureCounts_reverse/nofilter", out2 = "featureCounts_forward/nofilter"
    log: "logs/04-read_counts.log"
    conda: "alignment.yml"
    shell:
        "bash code/04-read_counts.sh {params.out1} {params.out2} {input.gff} {log} {input.bam}"

rule calculate_coverage:
    output:
        add_path_extension(all_input, "results/coverage", "perbase.cov")
    input:
        add_path_extension(all_input, "results/bam", "bam")
    conda: "circular.yml"
    shell:
        "bash code/05-coverage.sh"

rule picard_tools:
    output:
        pdf = add_path_extension(all_input, "results/picard", "pdf", "_insert_size_histogram"),
        txt = add_path_extension(all_input, "results/picard", "txt", "_insert_size_metrics")
    input:
        add_path_extension(all_input, "results/bam", "bam")
    conda: "alignment.yml"
    shell:
        "bash code/06-picard.sh"

rule get_TPM:
    output:
        per_sample = "results/TPM/TPM_per_sample.tsv",
        mean = "results/TPM/mean_TPM.tsv"
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts")
    params: os.getcwd()
    log: "logs/07-calculate_TPM.log"
    conda: "alignment.yml"
    script:
        "code/07-calculate_TPM.py"

rule filter_counts:
    output:
        counts = add_path_extension(all_input, "featureCounts_reverse/filtered", "featureCounts"),
        summary = add_path_extension(summary_input, "results/summary", "tsv", "_count_distribution")
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        gbff = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff" 
    conda: "alignment.yml"
    params: workdir = os.getcwd()
    log: "logs/08-prefilter_counts.py"
    script:
        "code/08-prefilter_counts.py"

rule parse_counts:
    output:
        avg_nofilter = "featureCounts_reverse/countfiles/nofilter/avg_gene_counts.tsv",
        avg_filtered = "featureCounts_reverse/countfiles/filtered/avg_gene_counts.tsv",
        meta = "featureCounts_reverse/countfiles/H3B1-04J_metadata.tsv",
        pos = "meta/gene_positions.tsv"
    input:
        TPM = "results/TPM/mean_TPM.tsv",
        meta = "meta/metadata.VF-3336.csv",
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        counts_filtered = add_path_extension(all_input, "featureCounts_reverse/filtered", "featureCounts")
    params: countdir = "featureCounts_reverse/nofilter", filtered_countdir = os.getcwd() + "/featureCounts_reverse/filtered"
    conda: "alignment.yml"
    log: "logs/09-filter_counts.log"
    script:
        "code/09-filter_counts.py"

# Add logging, update Python and bash variables so that they can be changed from this file
# and write rules for the saturation analysis and for the DE analysis to make sure that they generate the plots that I need

# R rules don't take the input from Snakemake, since I prefer not to install R with Conda
rule saturation:
    output:
        plots = expand("plots/saturation{extra}_k0.png", extra = ["_collapsed", ""])
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/filtered", "featureCounts"),
        meta = "featureCounts_reverse/countfiles/H3B1-04J_metadata.tsv"
    conda: "renv.yml"
    script:
        "code/10-saturation_analysis.R"

rule DE:   
    output:
        DE = expand("results/DE/{comparison}{ext}.csv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["_lfc1", ""]),
        PCA = expand("plots/pcaplot.{ext}", ext = ["png", "pdf"]),
        plots = expand("plots/{comparison}.{ext}", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["png", "ps"])
    input:
        counts = add_path_extension("featureCounts_reverse/filtered", "featureCounts")
    conda: "renv.yml"
    script:
        "code/11-dif_expression.R"

rule annotate_results:
    output: expand("results/DE/{comparison}{ext}_annotated.tsv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["", "_lfc1"])
    input:
        DE = expand("results/DE/{comparison}{ext}.csv", comparison = ["Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor", "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor"], ext = ["_lfc1", ""]),
        gbk = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff"
    conda: "alignment.yml"
    log: "logs/12-add_annotations.log"
    script: "12-add_annotations.py"
