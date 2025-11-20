infile=~/mucoid_project/adhesins/sequences/Muc_adhesins.faa
outbase=~/mucoid_project/adhesins/results/clustering/Muc_adhesins
mkdir -p $(dirname -- $outbase)
echo $infile
echo $outbase

cd-hit -i $infile -o $outbase -c 0.8 -aL 0.2 -aS 0.9 -g 1
