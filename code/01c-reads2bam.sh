index=$1
cores=$2
outputs=("${@:3:4}")
outdir=$(dirname -- ${outputs[0]})
inputs=("${@:7}")
count=0

echo "No. of threads: $cores"
echo "Index file: $index"
echo "Input files: ${inputs[@]}"
echo "Output BAM files: ${outputs[@]}"

echo "Creating output directory: $outdir"
mkdir -p $outdir

echo "Running BWA"
for file in ${inputs[@]};
do
outfile=$outdir/$(basename -- ${file%.fastq.gz}).bam

bwa mem -x pacbio -t $cores $index $file | samtools sort > $outfile
mv $outfile ${outputs[$count]}
echo 'Generated' ${outputs[$count]}
((count++))
done
