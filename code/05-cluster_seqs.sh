infile=~/adhesins/sequences/Muc_adhesins.faa
outbase=~/adhesins/results/clustering/Muc_adhesins
mkdir -p $(dirname -- $outbase)
echo $infile
echo $outbase

cd-hit -i $infile -o $outbase -c 0.9 -aL 0.2 -aS 0.9 -g 1
