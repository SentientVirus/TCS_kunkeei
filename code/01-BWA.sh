#!/bin/bash
####################################
###  RUN hisat2 + featureCounts  ###
###    (Samtools environment)    ###
####################################
workdir=$HOME'/snpseq00064'
path2files=$HOME'/snpseq00064/files/VF-3336/221006_M06455_0144_000000000-KMH8C'

index_file=$HOME'/Akunkeei_files/fna/H3B1-04J_genomic.fna'
gff_file=$HOME'/Akunkeei_files/gff/H3B1-04J_genomic.gff'

mkdir -p $workdir/trimmed_reads

# Genome indexing
cat $index_file > index/H3B1-04J.fna;
name='H3B1-04J'

hisat2-build $workdir/index/H3B1-04J.fna $workdir/index/H3B1-04J.index.hisat;
bwa index $workdir/index/H3B1-04J.fna > $workdir/index/H3B1-04J.bwt;

for sample in $path2files/Sample_VF-3336-H3B1-04J-*/*R1_001.fastq.gz;
do

# Read trimming
echo $sample, $sample ${sample::-15}R2_001.fastq.gz
echo $workdir/trimmed_reads/$(basename -- ${sample::-16})
skewer -m pe -Q 30 -t 12 -y TruSeq3-PE-2.fa $sample ${sample::-15}R2_001.fastq.gz -o $workdir/trimmed_reads/$(basename -- ${sample::-16})

# Alignment to generate bam files
bwa mem -t 12 $workdir/index/H3B1-04J.fna $workdir/trimmed_reads/$(basename -- ${sample::-16})-trimmed-pair1.fastq \
$workdir/trimmed_reads/$(basename -- ${sample::-16})-trimmed-pair2.fastq | samtools sort > $workdir/results/bam/$(basename -- ${sample::-16}).bam

# I checked that the experiment is reverse-stranded and kept only the code for the reverse strand
featureCounts -p -T 12 -M -C -s 2 -Q 10 -t gene -g locus_tag -a $gff_file \
-o $workdir/featureCounts_reverse/nofilter/$(basename -- ${sample::-16}).featureCounts $workdir/results/bam/$(basename -- ${sample::-16}).bam \
> $workdir/featureCounts_reverse/nofilter/$(basename -- ${sample::-16}).featureCounts

#featureCounts -p -T 12 -M -C -s 1 -Q 10 -t gene -g locus_tag -a index/$name.gff \
#-o featureCounts_forward/$(basename -- ${sample::-16}).log.vs.stat.featureCounts results/$(basename -- ${sample::-16}).bam \
#> featureCounts_forward/$(basename -- ${sample::-16}).log.vs.stat.featureCounts


done
