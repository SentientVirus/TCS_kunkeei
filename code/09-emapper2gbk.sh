#!/bin/bash

#############################################
###       Script to run emapper2gbk       ###
###     (annotations.yml environment)     ###
#############################################

fna=(${@:1:5})
faa=(${@:6:10})
gff=(${@:11:15})
annot=(${@:16:20})
output=(${@:21:25})
threads=${26}

echo $threads" threads will be used!"

for i in ${!fna[@]};
do
fna_file=${fna[$i]}
faa_file=${faa[$i]}
gff_file=${gff[$i]}
annot_file=${annot[$i]}
out_file=${output[$i]}
outdir=$(dirname -- $out_file)
mkdir -p $outdir

echo $i" Assembly file: "$fna_file
echo $i" Protein file: "$faa_file
echo $i" GFF file: "$gff_file
echo $i" EggNOG-mapper annotation: "$annot_file
echo $i" Output GenBank: "$out_file

emapper2gbk genomes -fn $fna_file -fp $faa_file -o $out_file -g $gff_file -gt CDS -n "Apilactobacillus kunkeei" -a $annot_file -c $threads --keep-gff-annotation 

done
