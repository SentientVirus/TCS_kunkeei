#!/bin/bash
##################################
###    RUN FastQC & MultiQC    ###
###    (read_QC environment)   ###
##################################

#Load inputs
inputs=${@:4}
cores=$1
outdir=$(dirname $2)
mqc_outdir=$(dirname $3)

mkdir -p $outdir #Create output directory if it doesn't exist

fastqc $inputs -o $outdir -t $cores #Run FastQC
multiqc -f $outdir/*.zip -o $mqc_outdir #Run MultiQC on FastQC output
