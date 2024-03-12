#!/bin/bash
##############################
###       RUN Skewer       ###
### (Samtools environment) ###
##############################

# Read variables from Snakemake
samples=${@:3}
outdir=$1
adapter=$2

mkdir -p $outdir

# Loop through samples
for sample in $samples;
do

# Print input and output filenames
echo $sample, ${sample/R1/"R2"}
echo $outdir/$(basename -- ${sample::-16})

# Trim reads
skewer -m pe -Q 30 -l 36 -t 12 -x $adapter -y $adapter $sample ${sample/R1/"R2"} -o $outdir/$(basename -- ${sample::-16})

done
