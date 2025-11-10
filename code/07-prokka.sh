#!/bin/bash

########################################
###  Script to annotate new genomes  ###
###  (Genome analysis environment)   ###
########################################

# Read variables from Snakemake
threads=$1
ref=$2
assemblies=${@:3:4}
outdir=$(dirname -- $7)

for assembly in $assemblies; #Loop through input assemblies
do

n=$(echo $assembly | cut -d"/" -f 3 | cut -d"." -f 1 | cut -d"i" -f 3) #Get isolate number
echo "Isolate number: "$n #Print isolate number

output_dir=${outdir::-3}/$n
echo "Output directory: "$output_dir

#Run Prokka
prokka $assembly --outdir $output_dir --prefix $n \
--force --addgenes --locustag AKUH3B104J --genus Apilactobacillus \
--species kunkeei --strain H3B1-04J --gram positive \
--usegenus --protein $ref --cpus $threads;

done
