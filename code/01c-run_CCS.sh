#!/bin/bash

#########################################
###         Script to run CCS         ###
### (genome_analysis.yml environment) ###
#########################################

infile=$1 #Input file
outfile=$2 #Output file

echo "Input: "$infile #Print path to input
echo "Output: "$outfile #Print path to output

echo "Run CCS"
ccs $infile $outfile
