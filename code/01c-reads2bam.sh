#!/bin/bash

########################################
###  Script to run BWA to align the  ###
###     reads back to the genome     ###
###      (Samtools environment)      ###
########################################

#Define input variables
index=$1
cores=$2
outputs=("${@:3:4}")
outdir=$(dirname -- ${outputs[0]})
inputs=("${@:7}")
count=0

#Print input variables
echo "No. of threads: $cores"
echo "Index file: $index"
echo "Input files: ${inputs[@]}"
echo "Output BAM files: ${outputs[@]}"
echo "Creating output directory: $outdir"

#Create output directory if it doesn't exist
mkdir -p $outdir

echo "Running BWA"
for file in ${inputs[@]}; #Loop through input files
do
outfile=$outdir/$(basename -- ${file%.fastq.gz}).bam #Define path to output file

#Run BWA and pipe the result to Samtools to convert it to BAM
bwa mem -x pacbio -t $cores $index $file | samtools sort > $outfile
mv $outfile ${outputs[$count]} #Move the output file to the desired path
echo 'Generated' ${outputs[$count]}
((count++)) #Increase the value of the variable to loop through outputs
done
