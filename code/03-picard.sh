#!/bin/bash
workdir=$HOME/snpseq00064
outdir=$workdir/results/picard
mkdir -p $outdir
for input in $workdir/results/bam/*.bam;
do
picard CollectInsertSizeMetrics \
I=$input \
O=$outdir/$(basename -- ${input%.bam})_insert_size_metrics.txt \
H=$outdir/$(basename -- ${input%.bam})_insert_size_histogram.pdf \
M=0.5
done
