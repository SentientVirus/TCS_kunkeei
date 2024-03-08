inputs=${@:4}
cores=$1
outdir=$(dirname $2)
mqc_outdir=$(dirname $3)

mkdir -p $outdir

fastqc $inputs -o $outdir -t $cores
multiqc -f $outdir -o $mqc_outdir
