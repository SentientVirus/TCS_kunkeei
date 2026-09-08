#!/bin/bash
##############################
###    RUN Picard Tools    ###
###    (alignment.yml)     ###
##############################

# Read inputs from Snakemake
outdir=$1
bam_files=${@:2}

# Create output directory if it doesn't exist
mkdir -p $outdir

# Loop through input files
for input in $bam_files;
do

# Run Picard Tools
picard CollectInsertSizeMetrics \
I=$input \
O=$outdir/$(basename -- ${input%.bam})_insert_size_metrics.txt \
H=$outdir/$(basename -- ${input%.bam})_insert_size_histogram.pdf \
M=0.5

done
