#!/bin/bash
############################################
###  Script to index genome with hisat2  ###
###        (Samtools environment)        ###
############################################

# Read variables from Snakemake
index_file=$1
outfile=$2
indir="$(dirname "${outfile}")"

# Create index directory if it does not exist
mkdir -p $indir

# Copy contents to the index file to the index directory
cat $index_file > $outfile;

# Index genome
bwa index $outfile
