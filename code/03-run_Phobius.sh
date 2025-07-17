#!/bin/bash
##############################
###       RUN Phobius      ###
###    (python_env.yml)    ###
##############################
infile=$1 #Input file with all protein sequences in strain H3B1-04J in FASTA format
outfile=$2 #Path to output file
logfile=$3 #Path to log file

mkdir -p $(dirname $outfile) #Create output directory if it doesn't exist
mkdir -p $(dirname $logfile) #Create log directory if it doesn't exist

phobius $infile > $outfile 2> $logfile #Run Phobius
