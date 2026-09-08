#!/bin/bash

#############################################
###  Script to run Circlator to set oriC  ###
###   (genome_analysis.yml environment)   ###
#############################################

# Read variables from Snakemake
inputs=("${@:1:4}")
outputs=("${@:5:4}")
outdir=$(dirname $5)
start_genes=$9
log=$10
> $log #Overwrite log file 

#Save paths to inputs and outputs to log file
echo "Inputs: $inputs" >> $log
echo "Outputs: $outputs" >> $log
echo "Output directory: $outdir" >> $log

# Create output directory if it doesn't exist
mkdir -p $outdir

# Run Circlator on each file
for i in {0..3};
do
input=${inputs[$i]} >> $log
output=${outputs[$i]} >> $log

echo "Processing $input, writing $output..." >> $log
circlator fixstart --genes_fa $start_genes $input ${output%.fasta} 2>> $log;
done
