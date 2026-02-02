#!/bin/bash

#########################################
### Script to demultiplex PacBio data ###
### (genome_analysis.yml environment) ###
#########################################

#Input files
bam=$1 #Multiplexed BAM file after CCS
barcode=$2 #FASTA file with barcode sequences
demux=${@:3:4} #Demultiplexed outputs (BAM format)
fastq=${@:7:4} #Demultiplexed outputs (FASTQ format)
outdir=$(dirname -- ${demux%% *}) #Output directory

echo "Create output directory if it doesn't exist: "$outdir
mkdir -p $outdir

echo "Multiplexed input: "$bam
echo "Path to file with barcodes: "$barcode
echo "Demultiplexed output (BAM): "$demux
echo "Demultiplexed output (fastq): "$fastq

#Convert output file lists to arrays
demux_arr=($demux)
fastq_arr=($fastq)

#Run Lima to demultiplex the CCS BAM file
lima $bam $barcode $outdir/$outdir.bam --same --split-bam

for i in "${!demux_arr[@]}"; do #Loop through output files by index
    echo "Converting "${demux_arr[$i]}" into "${fastq_arr[$i]}
    samtools bam2fq ${demux_arr[$i]} > ${fastq_arr[$i]} #Convert the files to FASTQ format
done
