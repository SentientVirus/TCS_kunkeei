#!/bin/bash

#############################################
###       Script to run emapper2gbk       ###
###     (annotations.yml environment)     ###
#############################################

#Load inputs from Snakemake
fna=(${@:1:5})
faa=(${@:6:10})
gff=(${@:11:15})
annot=(${@:16:20})
output=(${@:21:25})
threads=${26}

#Print number of threads
echo $threads" threads will be used!"

#Loop through inputs to generate outputs
for i in ${!fna[@]};
do
fna_file=${fna[$i]} #Retrieve all input files by index
faa_file=${faa[$i]}
gff_file=${gff[$i]}
annot_file=${annot[$i]}
out_file=${output[$i]}
outdir=$(dirname -- $out_file) #Retrieve the path to the output directory
mkdir -p $outdir #Make the output directory if it doesn't exist

#Print the paths to all inputs and outputs
echo $i" Assembly file: "$fna_file
echo $i" Protein file: "$faa_file
echo $i" GFF file: "$gff_file
echo $i" EggNOG-mapper annotation: "$annot_file
echo $i" Output GenBank: "$out_file

#Run emapper2gbk, -fn = genome fasta, -fp = protein fasta, -o output GenBank, -g = gff from Prokka, -gt = feature type, -n = organism name, -a = eggNOG annotation, -c = no. of threads
emapper2gbk genomes -fn $fna_file -fp $faa_file -o $out_file -g $gff_file -gt CDS -n "Apilactobacillus kunkeei" -a $annot_file -c $threads --keep-gff-annotation -go obo-db/go-basic.obo

done
