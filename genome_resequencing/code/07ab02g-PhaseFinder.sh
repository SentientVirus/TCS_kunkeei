#!/bin/bash

###################################
###  Script to run PhaseFinder  ###
###    (Samtools environment)   ###
###################################

# Load inputs from Snakemake
outdir=$(dirname $1)
genomes=${@:2}

# Create output directory if it doesn't exist
mkdir -p $outdir

# Loop through each genome and run PhaseFinder
for genome in $genomes;
do
prefix=$(basename "${genome}" | cut -d '.' -f1);
echo $genome
python PhaseFinder/PhaseFinder.py locate -f $genome -t $outdir/$prefix.tab -g 15 85 -p; 
done
