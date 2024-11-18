#!/bin/bash
infile=~/Akunkeei_files/faa/H3B1-04J_protein.faa #Input file with all protein sequences in strain H3B1-04J in FASTA format
outfile=~/proteomics/results/SignalP/H3B1-04J_SignalP.txt #Path to output file
logfile=~/proteomics/logs/signalp.log #Path to log file

mkdir -p $(dirname $outfile) #Create output directory if it doesn't exist
mkdir -p $(dirname $logfile) #Create log directory if it doesn't exist

~/signalp-4.1/signalp -f 'short' -t 'gram+' -v -l $logfile $infile > $outfile #Run SignalP
