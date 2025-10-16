#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 18:29:09 2025

Script to check the type of LPXTG in each gene.

@author: marina
"""

import os
import re
from Bio import SeqIO

workdir = os.path.expanduser('~') + '/adhesins'
indir = f'{workdir}/sequences/adhesins'
outdir = f'{workdir}/results/adhesin_lists'

adhesins = ['MucBP+LPXTG', 'MubB2+LPXTG']

infiles = [f'{indir}/{adhesin}.faa' for adhesin in adhesins]

for infile in infiles:
    type_dict = {}
    gene_type = os.path.basename(infile).replace('.faa', '')
    outfile = f'{outdir}/{gene_type}_LPXTG.tsv'
    with open(infile) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            LPXTG_motif = re.findall(r'LP(.{1,1})TG', str(record.seq))
            if len(LPXTG_motif) == 1:
                motif = LPXTG_motif[0]
            elif len(LPXTG_motif) == 0:
                motif = '-'
            else:
                motif = 'Several'
            type_dict[record.id] = motif
    with open(outfile, 'w') as out_handle:
        out_handle.write('Locus_tag\tLPXTG_type\n')
        {out_handle.write(f'{k}\t{v}\n') for (k, v) in type_dict.items()}