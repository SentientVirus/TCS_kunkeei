#configfile: "config.yml"

import os

comparisons = ["aggF_log_vs_aggS_log", "aggF_stat_vs_aggS_stat", "aggS_log_vs_aggS_stat", 
               "mucFD_stat_vs_mucF_log", "mucF_log_vs_mucS_log", "mucSD_stat_vs_mucS_log",
               "aggF_log_vs_mucF_log", "aggF_stat_vs_mucFD_stat", "aggS_log_vs_mucS_log", 
               "mucFD_stat_vs_mucF_stat", "mucF_log_vs_mucF_stat", "aggF_log_vs_aggF_stat",
               "aggF_stat_vs_mucF_stat", "aggS_stat_vs_mucSD_stat", "mucFD_stat_vs_mucSD_stat",
               "mucF_stat_vs_mucSD_stat"]

comparisons_DE = ["aggF_log_vs_aggS_log", "aggF_stat_vs_aggS_stat", "aggS_log_vs_aggS_stat",
               "mucFD_stat_vs_mucF_log", "aggF_log_vs_mucF_log", "aggF_stat_vs_mucFD_stat",
               "mucFD_stat_vs_mucF_stat", "mucF_log_vs_mucF_stat", "aggF_log_vs_aggF_stat",
               "aggF_stat_vs_mucF_stat", "aggS_stat_vs_mucSD_stat", "mucFD_stat_vs_mucSD_stat",
               "mucF_stat_vs_mucSD_stat"]


##Rule to generate the DAG
rule all:
    input:
        "results/TMH_predictions/TMH.tab",
        plots = expand("plots/{comparison}.png", comparison = comparisons_DE)

##Rule to parse the MS results
rule parse_MS:
    output:
        general = "files/parsed/all_comparisons.tsv",
        pairwise = expand("files/parsed/{comparison}.tsv", comparison = comparisons)
    input:
        exp1 = "files/MS-24-030_MaxQuant_results.xlsx",
        exp2 = "files/MS-24-038_MaxQuant_results.xlsx",
        gbk1 = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff",
        gbk2 = os.path.expanduser("~") + "/mucoid_project/ugc00027/results/annotations/emapper2gbk/reference.gbk"
    params: outdir = "files/parsed"
    conda: "pixi_proteome/default.yml"
    log: "logs/01-parse_MS.log"
    script:
        "code/01-parse_proteomics.py"

##Rule to add locus tag and SignalP information to the files (SignalP is run online)
rule add_SignalP:
    output:
        loci = expand("files/SP/{comparison}.tsv", comparison = comparisons),
        signalP = "results/SignalP/H3B1-04J_SignalP.tsv"
    input:
        infiles = expand("files/parsed/{comparison}.tsv", comparison = comparisons),
        signalP = "results/SignalP/prediction_results.txt" #Generated from the web server
    params: outdir = "files/loci"
    conda: "pixi_proteome/default.yml"
    log: "logs/02-add_SignalP.log"
    script:
        "code/02-add_SP.py"

##Rule to run DeepTMHMM
rule DeepTMHMM:
    output:
        "results/DeepTMHMM/deeptmhmm_results.md",
        "results/DeepTMHMM/predicted_topologies.3line",
        "results/DeepTMHMM/TMRs.gff3"
    input:
        proteins = os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa"
    params: outdir = "results/DeepTMHMM"
    conda: "pixi_proteome/default.yml"
    log: "logs/03-run_DeepTMHMM.log"
    script:
        "code/03-run_DeepTMHMM.py"

##Rule to run Phobius
rule Phobius:
    output:
        "results/Phobius/H3B1-04J_phobius.txt"
    input:
        os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa"
    conda: "pixi_proteome/default.yml"
    log: "logs/03-run_Phobius.log"
    shell:
        "bash code/03-run_Phobius.sh {input} {output} {log}"

##Rule to combine TMHMM predictions and incorporate them into the proteomics file
rule combine_TMHMM:
    output:
        common_pred = "results/TMH_predictions/TMH.tab",
        sample_pred = expand("files/SP/{comparison}_TMH.tsv", comparison = comparisons) 
    input:
        gbk = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff",
        Phobius = "results/Phobius/H3B1-04J_phobius.txt",
        DeepTMHMM = "results/DeepTMHMM/TMRs.gff3",
        sample_in = expand("files/loci/{comparison}.tsv", comparison = comparisons)
    conda: "pixi_proteome/default.yml"
    log: "logs/04-combine_TMH_predictions.log"
    script:
        "code/04-combine_TMH_predictions.py"


##Rule to generate Volcano plots in R
rule volcano:
    output:
        plots = expand("plots/{comparison}.png", comparison = comparisons_DE)
    input:
        infiles = expand("files/loci/{comparison}.tsv", comparison = comparisons_DE)
    conda: "pixi_proteome/renv.yml"
    log: "logs/05-volcano_plot.log"
    script:
        "code/05-volcano_proteomics.R"
