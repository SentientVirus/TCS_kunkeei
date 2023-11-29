#!/bin/bash
##############################
###   RUN featureCounts    ###
###    (alignment.yml)     ###
##############################

outdir1=$1
outdir2=$2

gff_file=$3

log=$4

samples=${@:5}

> $log

mkdir -p $outdir1
mkdir -p $outdir2

for sample in $samples;
do
echo $sample
echo $outdir1/$(basename -- ${sample::-4}).featureCounts

# Counts for the reverse strand
featureCounts -p -T 12 -M -C -s 2 -Q 10 -t gene -g locus_tag -a $gff_file \
-o $outdir1/$(basename -- ${sample::-4}).featureCounts $sample 2>> $log 3>> $log 

#Counts for the forward strand
featureCounts -p -T 12 -M -C -s 1 -Q 10 -t gene -g locus_tag -a $gff_file \
-o $outdir2/$(basename -- ${sample::-4}).featureCounts $sample 2>> $log 3>> $log

done
