indir="data/fixed_ori"
suffix="_genomic.fna"
outdir="results/PhaseFinder"
genomes=()

mkdir -p $outdir

for ((j=1; j<5; j++))
do
genomes+=("${indir}/sample${j}${suffix}")
done
genomes+=("../Akunkeei_files/fna/H3B1-04J${suffix}")

for ((i=0; i<${#genomes[@]}; i++))
do
prefix=$(basename "${genomes[i]}" | cut -d '.' -f1);
python PhaseFinder/PhaseFinder.py locate -f ${genomes[i]} -t $outdir/$prefix.tab -g 15 85 -p; 
done
