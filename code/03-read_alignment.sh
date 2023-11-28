#!/bin/bash
##################################
###          RUN BWA           ###
###  (Samtools environment)    ###
##################################


outdir=$1
samples=${@:4}
fna_file=$2
log=$3

mkdir -p $outdir

> $log

for sample in $samples;
do

echo "Processing: "$(basename -- ${sample::-20})  2>> $log

# Alignment to generate bam files
bwa mem -t 12 $fna_file $sample $(${sample}/pair1/"pair2") | samtools sort > $outdir/$(basename -- ${sample::-20}).bam 2>> $log 3>> $log

done
