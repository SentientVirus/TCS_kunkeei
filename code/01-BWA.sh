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

cat $index_file > index/H3B1-04J.fna;
name='H3B1-04J'

hisat2-build $workdir/index/H3B1-04J.fna $workdir/index/H3B1-04J.index.hisat;
bwa index $workdir/index/H3B1-04J.fna > $workdir/index/H3B1-04J.bwt;

#for strain in $path2files/Sample_*;
#do
#name=$(echo $strain | grep -oP '(?<=4-).*');

#file=$(echo /home/marina/MolevoUU/Bronze/molevol/LAB/annotations_220914/fna/$name.fna);
#cp $file index/H3B1-04J.fna;
#hisat2-build index/$name.fna index/$name.index.hisat2;

for sample in $path2files/Sample_VF-3336-H3B1-04J-*/*R1_001.fastq.gz;
do
echo $sample, $sample ${sample::-15}R2_001.fastq.gz
echo $workdir/trimmed_reads/$(basename -- ${sample::-16})
skewer -m pe -Q 30 -t 12 -y TruSeq3-PE-2.fa $sample ${sample::-15}R2_001.fastq.gz -o $workdir/trimmed_reads/$(basename -- ${sample::-16})

#bwa index index/$name.fna > index/$name.bwt
bwa mem -t 12 $workdir/index/H3B1-04J.fna $workdir/trimmed_reads/$(basename -- ${sample::-16})-trimmed-pair1.fastq \
$workdir/trimmed_reads/$(basename -- ${sample::-16})-trimmed-pair2.fastq | samtools sort > $workdir/results/bam/$(basename -- ${sample::-16}).bam

#This is the right strand
featureCounts -p -T 12 -M -C -s 2 -Q 10 -t gene -g locus_tag -a $gff_file \
-o $workdir/featureCounts_reverse/nofilter/$(basename -- ${sample::-16}).featureCounts $workdir/results/bam/$(basename -- ${sample::-16}).bam \
> $workdir/featureCounts_reverse/nofilter/$(basename -- ${sample::-16}).featureCounts

#featureCounts -p -T 12 -M -C -s 1 -Q 10 -t gene -g locus_tag -a index/$name.gff \
#-o featureCounts_forward/$(basename -- ${sample::-16}).log.vs.stat.featureCounts results/$(basename -- ${sample::-16}).bam \
#> featureCounts_forward/$(basename -- ${sample::-16}).log.vs.stat.featureCounts


done

#Add stranded option (-s 2) and filter low counts (after counting) and check average fragment length (76bp, same for all sequences)
