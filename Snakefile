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
