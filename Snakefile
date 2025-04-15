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
        exp2 = "files/MS-24-038_MaxQuant_results.xlsx"
    params: outdir = "files/parsed"
    conda: "envs/python_env.yml"
    log: "logs/01-index_genome.log"
    shell:
        "bash code/01-index_genome.sh {input} {output.fna} 2> {log}"
