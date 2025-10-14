indir=~/adhesins/sequences/adhesins
outdir=~/adhesins/alignments/adhesins
logfile=~/adhesins/logs/04-align_sequences.log
treedir=~/adhesins/trees/adhesins
threads=8

mkdir -p $(dirname -- $logfile)
mkdir -p $outdir
mkdir -p $treedir

> $logfile

for file in $indir/*.faa;
do
outfile=$(basename -- $file)
outfile=$outdir/$(echo $outfile | cut -d'.' -f 1).mafft.$(echo $file | cut -d'.' -f 2)

echo 'Input: '$file
echo 'Output: '$outfile
echo 'Log: '$logfile

mafft-linsi --thread $threads $file > $outfile 2>> $logfile;
iqtree -nt AUTO -ntmax $threads -s $outfile -st AA -msub nuclear -bb 1000 -bnni >> $logfile
mv $outfile.* $treedir
done
