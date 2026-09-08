#!/bin/bash
##################################
###          RUN BWA           ###
###  (Samtools environment)    ###
##################################

# Read inputs from Snakemake
outdir=$1
samples=${@:4}
fna_file=$2
log=$3

# Create output directory if it doesn't exist
mkdir -p $outdir

# Create or empty log file
> $log

# Loop through input files
for sample in $samples;
do

# Print name of the file to log
echo "Processing: "$(basename -- ${sample::-20}) >> $log 2>> $log

# Alignment to generate bam files
bwa mem -t 12 $fna_file $sample ${sample/pair1/"pair2"} 2>> $log | samtools sort > $outdir/$(basename -- ${sample::-20}).bam 2>&1 2>> $log

done
