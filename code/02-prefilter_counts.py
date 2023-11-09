#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:32:56 2023

This is a script to filter counts that removes ribosomal genes

@author: marina
"""
# =============================================================================
# 0. Importing packages
# =============================================================================
from os.path import expanduser, exists, isfile, basename
from os import makedirs, listdir
from Bio import GenBank
import pandas as pd

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================
home = expanduser('~')

workdir = f'{home}/snpseq00064'

count_dir = f'{workdir}/featureCounts_reverse'
outdir = f'{count_dir}/filtered'
metadir = f'{workdir}/results/summary'

if not exists(outdir):
    makedirs(outdir)
    
if not exists(metadir):
    makedirs(metadir)
    
gbk_file = f'{home}/Akunkeei_files/gbff/H3B1-04J_genomic.gbff'

# =============================================================================
# 2. Read GenBank file to get locus tags of genes annotated as ribosomal RNAs
# or tRNAs
# =============================================================================
loctags_by_type = {'rRNA': [], 'tRNA': [], 'misc_RNA': [], 'transposon': []}
RNA_loctags = []

with open(gbk_file) as handle:
    for record in GenBank.parse(handle):
        #print(record.accession)
        #exclude = ['RNA', 'transposon', 'transposase']
        for feature in record.features:
            if 'RNA' in feature.key or (len(feature.qualifiers) > 5 and 'transpos' in feature.qualifiers[5].value): #or (len(feature.qualifiers) > 1 and 'RNA' in feature.qualifiers[1].value)
                loctag = feature.qualifiers[0].value.replace('"', '')
                print(f'{loctag} in chromosome {record.accession[0]} appended to list.', end = ' ')
                RNA_loctags.append(loctag)
                if 'rRNA' in feature.key:
                    type_seq = 'rRNA'
                elif 'tRNA' in feature.key or 'ssrA' in loctag:
                    type_seq = 'tRNA'
                elif 'RNA' in feature.key or 'RNA' in feature.qualifiers[1].value:
                    type_seq = 'misc_RNA'
                else: type_seq = 'transposon'
                print(f'Element is a {type_seq}.')
                loctags_by_type[type_seq].append(loctag)
                
        
# =============================================================================
# 3. Read count file and filter out the RNA genes, create 
# =============================================================================
type_count = {}
for file in listdir(count_dir):
    file_path = f'{count_dir}/{file}'
    if file.endswith('featureCounts'):
        with open(f'{file_path}') as counts:
            df = pd.read_csv(counts, sep = '\t', header = 1, index_col = 0)
            total_reads = df.iloc[:, -1:].sum()[0]
            for seqtype in ['rRNA', 'tRNA', 'misc_RNA', 'transposon']:
                type_count[seqtype] = df[df.index.isin(loctags_by_type[seqtype])].iloc[:, -1:].sum()[0]
            df = df[~df.index.isin(RNA_loctags)]
            # df.columns = df.iloc[0]
        isolate = f'I{file[17:19]}_{file[20:21]}_{file[24:27].replace("_", "")}'
        with open(f'{outdir}/{file}', 'w') as outfile:
            df.to_csv(outfile, sep = '\t')
        with open(f'{metadir}/{isolate}_count_distribution.tsv', 'w') as metadat:
            metadat.write('feature_type\tcount\n')
            [metadat.write(f'{k}\t{v}\n') for k, v in type_count.items()]
            