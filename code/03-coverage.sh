#!/bin/bash
workdir=$HOME/snpseq00064
outdir=$workdir/results/coverage

mkdir -p $outdir
for bam in $workdir/results/bam/*.bam
do
    echo 'Processing '$(basename -- $bam)
    out=$workdir'/results/coverage/'$(basename -- ${bam%.bam})'.perbase.cov'
    genomeCoverageBed -ibam $bam -d > $out
done
