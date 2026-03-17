#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 17:03:59 2026

Script to merge the InterProScan results for each file into a single file.
The script is numbered as "03" to remind the user that there is an intermediate
step, which is running InterProScan.

@author: Marina Mota-Merlo
"""
import os
import pandas as pd

workdir = os.path.expanduser('~') + '/mucoid_project/snpseq00064'

indir = f'{workdir}/results/InterProScan'

col_list = ['Protein accession', 'Sequence MD5 digest', 'Sequence length',
              'Analysis', 'Signature accession', 'Signature description',
              'Start location', 'Stop location', 'Score', 'Status', 'Date',
              'InterPro accession', 'InterPro description', 'GO annotations',
              'Pathways annotations']

infiles = [f'{indir}/{file}' for file in os.listdir(indir) if 'proteins' in file and file.endswith('.tsv') and 'complete' not in file]
infiles = sorted(infiles)
outfile = f'{indir}/H3B1-04J_complete_proteins.tsv'
df_list = []

for file in infiles:
    df = pd.read_csv(file, sep = '\t', header = None)
    df_list.append(df)
    
df = pd.concat(df_list, ignore_index = True, sort = True, axis = 0)

df.columns = col_list

df[['Protein accession', 'Locus tag']] = df['Protein accession'].str.split('/', n = 1, expand = True)

df_final = df[[col_list[0], 'Locus tag'] + col_list[1:]]

df_final.sort_values(df_final.columns[1], axis = 0, inplace = True, ignore_index = True)

with open(outfile, 'w') as handle:
    df_final.to_csv(handle, sep ='\t', index = False)