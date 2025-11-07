#Script to generate a DAG plot
snakemake --dag | grep -v "\-> 0\|0\[label = \"all\"" | dot -Tpdf -o dag.pdf
