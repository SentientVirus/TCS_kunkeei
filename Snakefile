#configfile: "config.yml"

import os

comparisons = ["inhF_log_vs_inhS_log", "inhF_stat_vs_inhS_stat", "inhS_log_vs_inhS_stat", 
               "mucFD_stat_vs_mucF_log", "mucF_log_vs_mucS_log", "mucS_log_vs_mucSD_stat",
               "inhF_log_vs_mucF_log", "inhF_stat_vs_mucFD_stat", "inhS_log_vs_mucS_log", 
               "mucFD_stat_vs_mucF_stat", "mucF_stat_vs_mucF_log", "inhF_stat_vs_inhF_log",
               "inhF_stat_vs_mucF_stat", "inhS_stat_vs_mucSD_stat", "mucFD_stat_vs_mucSD_stat",
               "mucF_stat_vs_mucSD_stat"]

##Rule to index the reference genome of strain H3B1-04J
rule parse_MS:
    output:
        expand("files/parsed/{comparison}.tsv", comparison = comparisons)
    input:
        exp1 = "files/MS-24-030_MaxQuant_results.xlsx",
        exp2 = "files/MS-24-038_MaxQuant_results.xlsx",
        gbk = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff"
    params: outdir = "files/parsed"
    conda: "envs/python_env.yml"
    log: "logs/01-parse_MS.log"
    script:
        "code/01-parse_proteomics.py"

##Rule to run SignalP
rule run_SignalP:
    output:
        "results/SignalP/H3B1-04J_SignalP.txt"
    input: 
        os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa"
    conda: "envs/python_env.yml"
    log: "logs/01-run_SignalP.log"
    shell:
        "code/01-run_SignalP.sh {input} {output} {log}"

##Rule to add locus tag and SignalP information to the files
rule add_SignalP:
    output:
        loci = expand("files/loci/{comparison}.tsv", comparison = comparisons),
        signalP = "results/SignalP/H3B1-04J_SignalP.tsv"
    input:
        infiles = expand("files/parsed/{comparison}.tsv", comparison = comparisons),
        gbk = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff",
        signalP = "results/SignalP/H3B1-04J_SignalP.txt"
    params: outdir = "files/loci"
    conda: "envs/python_env.yml"
    log: "logs/02-add_SignalP.log"
    script:
        "code/02-add_SP_loctags.py"

##Rule to run DeepTMHMM
rule DeepTMHMM:
    output:
        "results/DeepTMHMM/deeptmhmm_results.md",
        "results/DeepTMHMM/predicted_topologies.3line",
        "results/DeepTMHMM/TMRs.gff3"
    input:
        proteins = os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa"
    params: outdir = "results/DeepTMHMM"
    conda: "envs/python_env.yml" #OBS! Install and activate the environment
    log: "logs/03-run_DeepTMHMM.log"
    script:
        "code/03-run_DeepTMHMM.py"

##Rule to run Phobius
rule Phobius:
    output:
        "results/Phobius/H3B1-04J_phobius.txt"
    input:
        os.path.expanduser("~") + "/Akunkeei_files/faa/H3B1-04J_protein.faa"
    conda: "envs/python_env.yml"
    log: "logs/03-run_Phobius.log"
    shell:
        "bash code/03-run_Phobius.sh {input} {output} {log}"

##Rule to combine TMHMM predictions and incorporate them into the proteomics file
rule combine_TMHMM:
    output:
        common_pred = "results/TMH_predictions/TMH.tab",
        sample_pred = expand("files/loci/{comparison}_TMH.tsv", comparison = comparisons) 
    input:
        gbk = os.path.expanduser("~") + "/Akunkeei_files/gbff/H3B1-04J_genomic.gbff",
        Phobius = "results/Phobius/H3B1-04J_phobius.txt",
        DeepTMHMM = "results/DeepTMHMM/TMRs.gff3",
        sample_in = expand("files/loci/{comparison}.tsv", comparison = comparisons)
    conda: "envs/python_env.yml"
    log: "logs/04-combine_TMH_predictions.log"
    script:
        "code/04-combine_TMH_predictions.py"
