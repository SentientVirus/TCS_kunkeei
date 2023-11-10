configfile: "config.yml"
import os

isol_code = {"01-F-1": "S5", "01-F-2": "S13", "01-F-3": "S21", "01-F-4": "S29", "01-F-5": "S37", "01-S-1": "S1", "01-S-2": "S9", 
    "01-S-3": "S17", "01-S-4": "S25", "01-S-5": "S33", "02-F-1": "S6", "02-F-2": "S14", "02-F-3": "S22", "02-F-4": "S30", "02-F-5": "S38",
    "02-S-1": "S2", "02-S-2": "S10", "02-S-3": "S18", "02-S-4": "S26", "02-S-5": "S34", "09-F-1": "S7", "09-F-2": "S15", "09-F-3": "S23", 
    "09-F-4": "S31", "09-F-5": "S39", "09-S-1": "S3", "09-S-2": "S11", "09-S-3": "S19", "09-S-4": "S27", "09-S-5": "S35", "10-F-1": "S8",
    "10-F-2": "S16", "10-F-3": "S24", "10-F-4": "S32", "10-F-5": "S40"}

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

rule align2fna:
    output:
        #trimmed_reads = add_path_extension(),
        outbam = add_path_extension(all_input, "results/bam", "bam"),
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts") 
    input:
        input_def(all_input, path = "files/VF-3336/221006_M06455_0144_000000000-KMH8C", extension = "fastq.gz")
    conda: "alignment.yml"
    shell:
        "bash code/01-BWA.sh"

rule calculate_coverage:
    output:
        add_path_extension(all_input, "results/coverage", "perbase.cov")
    input:
        add_path_extension(all_input, "results/bam", "bam")
    conda: "circular.yml"
    shell:
        "bash code/02-coverage.sh"

rule picard_tools:
    output:
        pdf = add_path_extension(all_input, "results/picard", "pdf", "_insert_size_histogram"),
        txt = add_path_extension(all_input, "results/picard", "txt", "_insert_size_metrics")
    input:
        add_path_extension(all_input, "results/bam", "bam")
    conda: "alignment.yml"
    shell:
        "bash code/03-picard.sh"

#rule get_TPM:
#    output:
#        ent = add_path_extension(all_input, "results/TPM", "ent", "_genes"),
#        out = add_path_extension(all_input, "results/TPM", "out", "_genes"),
#        uni = add_path_extension(all_input, "results/TPM", "uni", "_genes")
#    input:
#        bam = add_path_extension(all_input, "results/bam", "bam"),
#        gtf = "H3B1-04J.gtf"
#    conda: "circular.yml"
#    shell:
#        "bash code/04-TPM.sh"

rule get_TPM:
    output:
        per_sample = "results/TPM/TPM_per_sample.tsv",
        mean = "results/TPM/mean_TPM.tsv"
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts")
    params: os.getcwd()
    conda: "alignment.yml"
    script:
        "code/04-calculate_TPM.py"

rule filter_counts:
    output:
        counts = add_path_extension(all_input, "featureCounts_reverse/filtered", "featureCounts"),
        summary = add_path_extension(summary_input, "results/summary", "tsv", "_count_distribution")
    input:
        counts = add_path_extension(all_input, "featureCounts_reverse/nofilter", "featureCounts"),
        gbff = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff" 
    conda: "alignment.yml"
    params: workdir = os.getcwd()
    script:
        "code/05-prefilter_counts.py"

# Maybe here create a filtered version as well, without RNAs
#rule parse_TPM:
#    output:
#        "results/TPM/mean_TPM.tsv"
#    input:
#        ent = add_path_extension(all_input, "results/TPM", "ent", "_genes"),
#        out = add_path_extension(all_input, "results/TPM", "out", "_genes"),
#        uni = add_path_extension(all_input, "results/TPM", "uni", "_genes")
#    conda: "alignment.yml"
#    script:
#        "code/06-parse_TPM.py"

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
    params: countdir = os.getcwd() + "/featureCounts_reverse/nofilter", filtered_countdir = os.getcwd() + "/featureCounts_reverse/filtered"
    conda: "alignment.yml"
    script:
        "code/06-filter_counts.py"

# Add logging, update Python and bash variables so that they can be changed from this file
# and write rules for the saturation analysis and for the DE analysis to make sure that they generate the plots that I need
