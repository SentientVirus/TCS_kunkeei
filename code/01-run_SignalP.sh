#!/bin/bash
infile=~/Akunkeei_files/faa/H3B1-04J_protein.faa
outfile=results/SignalP/H3B1-04J_SignalP.txt

mkdir -p $(dirname $outfile)

../signalp-4.1/signalp -f 'short' -t 'gram+' -v -l logs/signalp.log $infile > $outfile
