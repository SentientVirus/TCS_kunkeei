outdir=$(dirname $1)
genomes=${@:2}

mkdir -p $outdir

for genome in $genomes;
do
prefix=$(basename "${genome}" | cut -d '.' -f1);
echo $genome
python PhaseFinder/PhaseFinder.py locate -f $genome -t $outdir/$prefix.tab -g 15 85 -p; 
done
