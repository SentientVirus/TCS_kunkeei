#!/bin/bash
workdir=$HOME/snpseq00064
for input in $workdir/results/bam/*.bam;
do
picard CollectInsertSizeMetrics \
I=$input \
O=results/picard/$(basename -- ${input%.bam})_insert_size_metrics.txt \
H=results/picard/$(basename -- ${input%.bam})_insert_size_histogram.pdf \
M=0.5
done
