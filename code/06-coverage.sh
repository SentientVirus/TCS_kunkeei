#!/bin/bash
##############################
###   RUN genomeCoverage   ###
###     (circular.yml)     ###
##############################

# Read inputs from Snakemake
outdir=$1
bam_files=${@:2}

# Create output directory if it doesn't exist
mkdir -p $outdir

# Loop through input files
for bam in $bam_files;
do

# Print input file name to log
echo 'Processing '$(basename -- $bam) >&2

# Define output
out=$outdir/$(basename -- ${bam%.bam})'.perbase.cov'

# Calculate coverage
genomeCoverageBed -ibam $bam -d > $out

done
