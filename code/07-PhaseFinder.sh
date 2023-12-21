outdir=$(dirname $1)
genomes=${@:2}

mkdir -p $outdir

for ((i=0; i<${#genomes[@]}; i++))
do
prefix=$(basename "${genomes[i]}" | cut -d '.' -f1);
python PhaseFinder/PhaseFinder.py locate -f ${genomes[i]} -t $outdir/$prefix.tab -g 15 85 -p; 
done
