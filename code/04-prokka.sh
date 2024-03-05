#!/bin/bash

########################################
###  Script to annotate new genomes  ###
###  (Genome analysis environment)   ###
########################################

# Read variables from Snakemake
threads=$1
ref=$2
assemblies=${@:3}

for assembly in $assemblies;
do
n=$(echo $assembly | cut -d"/" -f 3 | cut -d"." -f 1 | cut -d"i" -f 3)
echo $n
prokka $assembly --outdir results/annotations --prefix sample$n \
--force --addgenes --genus Apilactobacillus \
--species kunkeei --strain H3B1-04J --gram positive \
--usegenus --protein $ref --cpus $threads;
done
