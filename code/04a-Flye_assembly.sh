#!/bin/bash

###############################################
### Script to generate assemblies with Flye ###
###        (assembly.yml environment)       ###
###############################################

infiles=${@:1:4} #Input files
threads=$9 #No. of threads
echo "Inputs: "$infiles #Print input files

echo "Loop through inputs to generate assemblies"
i=0 #Index for the loop
for file in $infiles; #Loop through input files
do
echo "Index: "$i #Print the loop index
outfile=${@:5+$i:1} #Retrieve the output file
outdir=$(dirname -- $outfile) #Retrieve the output directory
echo "Outdir: "$outdir #Print output directory
echo "File: "$file #Print output file
i=$(expr $i + 1) #Increase the index
mkdir -p $outdir #Create the output directory if it does not exist
flye --pacbio-raw $file --genome-size 1.5m --threads $threads --out-dir $outdir --meta #Run Flye
done
