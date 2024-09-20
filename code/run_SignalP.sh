#!/bin/bash
infile=~/Akunkeei_files/faa/H3B1-04J_protein.faa
outfile=results/H3B1-04J_SignalP.tsv
outplots=plots/SignalP
mkdir -p $(dirname $outfile)
mkdir -p $outplots
../signalp-4.1/signalp -f 'short' -g 'png+eps' -t 'gram+' -v -l logs/signalp.log $infile > $outfile
mv *.png $outplots
mv *.eps $outplots
