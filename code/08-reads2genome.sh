outdir=$HOME/ugc00027/bam_files
indir=$HOME/snpseq00064
mkdir -p $outdir
for file in $HOME/ugc00027/subreads/*/*;
do
bwa mem -x pacbio -t 24 $indir/index/H3B1-04J.fna $file | samtools sort > $outdir/$(basename -- ${file%.fastq.gz}).bam
done

mv $outdir/*1002* 001.bam
mv $outdir/*1003* 002.bam
mv $outdir/*1008* 009.bam
mv $outdir/*1010* 010.bam
