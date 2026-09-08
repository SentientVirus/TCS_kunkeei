#!/bin/bash
##############################
###   RUN featureCounts    ###
###    (alignment.yml)     ###
##############################

# Read inputs from Snakemake
outdir1=$1
outdir2=$2
gff_file=$3
log=$4

samples=${@:5}

# Create or empty log file
> $log

# Create output directories if they don't exist
mkdir -p $outdir1
mkdir -p $outdir2

# Loop through samples
for sample in $samples;
do

# Print input and output filenames
echo $sample
echo $outdir1/$(basename -- ${sample::-4}).featureCounts

# Counts for the reverse strand
featureCounts -p -T 12 -M -C -s 2 -Q 10 -t gene -g locus_tag -a $gff_file \
-o $outdir1/$(basename -- ${sample::-4}).featureCounts $sample 2>> $log 3>> $log 

# Counts for the forward strand
featureCounts -p -T 12 -M -C -s 1 -Q 10 -t gene -g locus_tag -a $gff_file \
-o $outdir2/$(basename -- ${sample::-4}).featureCounts $sample 2>> $log 3>> $log

done
