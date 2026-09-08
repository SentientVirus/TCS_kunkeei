faa=~/Akunkeei_files/faa
outpath=~/mucoid_project/adhesins/interproscan
logpath=~/mucoid_project/adhesins/logs/01-interproscan

mkdir -p $outpath
mkdir -p $logpath

for file in $faa/*protein.faa;
do

outbase=$(basename -- ${file%_*})
outfile=$outpath/$outbase.tsv
logfile=$logpath/$outbase.log

echo Input: $file > $logfile
echo Output: $outfile >> $logfile
echo Log: $logfile >> $logfile

echo Starting run... >> $logfile
~/interproscan-5.59-91.0/interproscan.sh -i $file -f tsv -o $outfile -dp -cpu 12 2>> $logfile >> $logfile;

done

