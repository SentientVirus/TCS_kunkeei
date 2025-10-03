infile=~/adhesins/sequences/Muc_adhesins.faa
infile2=~/adhesins/sequences/Muc_adhesins_repset.faa
infiles=($infile $infile2)
logfile=~/adhesins/logs/04-align_sequences.log
treedir=~/adhesins/trees
threads=8

mkdir -p $(dirname -- $logfile)
mkdir -p $treedir

> $logfile

for file in ${infiles[@]};
do
outfile=$(echo $file | cut -d'.' -f 1).mafft.$(echo $file | cut -d'.' -f 2)

echo 'Input: '$infile
echo 'Output: '$outfile
echo 'Log: '$logfile

mafft-linsi --thread $threads $file > $outfile 2>> $logfile;
iqtree -nt AUTO -ntmax $threads -s $outfile -st AA -msub nuclear -bb 1000 -bnni >> $logfile
mv $outfile.* $treedir
done
