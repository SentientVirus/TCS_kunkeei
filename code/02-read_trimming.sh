#!/bin/bash
##############################
###       RUN Skewer       ###
### (Samtools environment) ###
##############################
path2files=$HOME'/snpseq00064/files/VF-3336/221006_M06455_0144_000000000-KMH8C'
samples=${@:3}
outdir=$1
adapter=$2

#echo $samples
#echo $outdir
#echo $adapter

mkdir -p $outdir

# Loop through samples
for sample in $samples;
do

# Print input and output filenames
echo $sample, ${sample/R1/"R2"}
echo $outdir/$(basename -- ${sample::-16})

# Trim reads
skewer -m pe -Q 30 -t 12 -y $adapter $sample ${sample/R1/"R2"} -o $outdir/$(basename -- ${sample::-16})

done
