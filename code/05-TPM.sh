workdir=$HOME/snpseq00064
outdir=$workdir/results/TPM
indir=$workdir/results/bam

mkdir -p $outdir 
cd $outdir

TPMCalculator -g $workdir/H3B1-04J.gtf -d $indir -c 76
