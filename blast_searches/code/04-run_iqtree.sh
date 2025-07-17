infile=$HOME/mucoid_project/blast_searches/fasta/blastp_formatted.mafft.faa
logfile=$HOME/mucoid_project/blast_searches/logs/iqtree.log

mkdir -p $(dirname -- $logfile)
> $logfile

iqtree -nt AUTO -ntmax 24 -s $infile -st AA -msub nuclear -bb 1000 -bnni >> $logfile
