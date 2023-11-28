#!/bin/bash
##############################
###       RUN Skewer       ###
### (Samtools environment) ###
##############################

samples=${@:4}
outdir=$1
adapter=$2
log=$3

echo $samples
echo $outdir
echo $adapter

mkdir -p $outdir

> $log

# Loop through samples
for sample in $samples;
do

# Print input and output filenames
echo $sample, $sample ${sample::-15}R2_001.fastq.gz
echo $outdir/$(basename -- ${sample::-16})

# Trim reads
skewer -m pe -Q 30 -t 12 -y $adapter $sample ${sample::-15}R2_001.fastq.gz -o $outdir/$(basename -- ${sample::-16}) 2>> $log 3>>$log

done
