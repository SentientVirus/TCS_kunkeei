#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 17:52:07 2026

@author: marina
"""

import os, re
import subprocess
from Bio import SeqIO

workdir = os.path.expanduser('~') + '/mucoid_project/ugc00027/blast_searches'
log = f'{workdir}/logs/aln_indentity.log'
infile = f'{workdir}/ref_aln/RR-TF.faa'
out_mafft = infile.replace('.faa', '.mafft.faa')
out_needle = infile.replace('.faa', '_needle.tsv')
sep_dir = f'{workdir}/ref_aln/separated_seqs'
in_prots = [f'{sep_dir}/{file}' for file in os.listdir(sep_dir) if file.endswith('.faa')]
t = 4
p = re.compile('\((.*)\)')
q = re.compile('Score: (.*)')

def run_mafft(infile, outfile, cores):
    command = f'mafft-linsi --thread {cores} {infile} > {outfile} 2> {log}'
    subprocess.run(command, shell = True)

run_mafft(infile, out_mafft, t)

with open(infile) as handle:
    rec_desc = {record.id: record.description for record in SeqIO.parse(handle, 'fasta')}

with open(out_needle, 'w') as out:
    out.write('seq1\tseq2\t%identity\t%similarity\t%gaps\tscore\n')
    

for i in range(len(in_prots)-1):
    aseq_name = os.path.basename(in_prots[i]).replace('.faa', '')
    for j in range(i + 1, len(in_prots)):
        bseq_name = os.path.basename(in_prots[j]).replace('.faa', '')
        outfile = f'{sep_dir}/{aseq_name}_vs_{bseq_name}.txt'
        command = f'needle -aseq {in_prots[i]} -bseq {in_prots[j]} -outfile stdout -sprotein 1 -aformat simple -gapopen 10 -gapextend 0.5 2>> {log}'
        result = subprocess.run(command, shell = True, capture_output = True)
        stdout_var = str(result.stdout).replace('\\n', '\n').split('\n')
        ident = p.search(stdout_var[25]).group(1).replace('%', '')
        simil = p.search(stdout_var[26]).group(1).replace('%', '')
        gaps = p.search(stdout_var[27]).group(1).replace('%', '')
        score = q.search(stdout_var[28]).group(1)
        with open(out_needle, 'a') as out:
            out.write(f'{aseq_name}\t{bseq_name}\t{ident}\t{simil}\t{gaps}\t{score}\n')
        print(aseq_name, bseq_name)
        print(ident, simil, gaps, score)


    