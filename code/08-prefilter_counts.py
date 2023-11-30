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

import os
import logging, traceback
from Bio import GenBank
import pandas as pd

# =============================================================================
# 0. Logging
# =============================================================================

logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================

workdir = snakemake.params.workdir

count_dir = os.path.dirname(snakemake.input.counts[0])
outdir = os.path.dirname(snakemake.output.counts[0])
metadir = os.path.dirname(snakemake.output.summary[0])

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
if not os.path.exists(metadir):
    os.makedirs(metadir)
    
gbk_file = snakemake.input.gbff

# =============================================================================
# 2. Read GenBank file to get locus tags of genes annotated as ribosomal RNAs
# or tRNAs
# =============================================================================
loctags_by_type = {'rRNA': [], 'tRNA': [], 'misc_RNA': [], 'transposon': []}
RNA_loctags = []

with open(gbk_file) as handle:
    for record in GenBank.parse(handle):
        for feature in record.features:
            if 'RNA' in feature.key or (len(feature.qualifiers) > 5 and 'transpos' in feature.qualifiers[5].value): #or (len(feature.qualifiers) > 1 and 'RNA' in feature.qualifiers[1].value)
                loctag = feature.qualifiers[0].value.replace('"', '')
                check = True
                if 'ssrA' in loctag:
                    loctag = feature.qualifiers[1].value.replace('"', '')
                    type_seq = 'tRNA'
                    check = False
                print(f'{loctag} in chromosome {record.accession[0]} appended to list.', end = ' ')
                RNA_loctags.append(loctag)
                if check and 'rRNA' in feature.key:
                    type_seq = 'rRNA'
                elif check and 'tRNA' in feature.key:
                    type_seq = 'tRNA'
                elif check and ('RNA' in feature.key or 'RNA' in feature.qualifiers[1].value):
                    type_seq = 'misc_RNA'
                elif check:
                    type_seq = 'transposon'
                print(f'Element is a {type_seq}.')
                loctags_by_type[type_seq].append(loctag)
                
        
# =============================================================================
# 3. Read count file and filter out the RNA genes, create filtered count files
# and files with the data of RNA genes 
# =============================================================================
type_count = {}
for file in os.listdir(count_dir):
    file_path = f'{count_dir}/{file}'
    if file.endswith('featureCounts'):
        with open(f'{file_path}') as counts:
            df = pd.read_csv(counts, sep = '\t', header = 1, index_col = 0)
            total_reads = df.iloc[:, -1:].sum()[0]
            for seqtype in ['rRNA', 'tRNA', 'misc_RNA', 'transposon']:
                type_count[seqtype] = df[df.index.isin(loctags_by_type[seqtype])].iloc[:, -1:].sum()[0]
            df = df[~df.index.isin(RNA_loctags)]
        isolate = f'I{file[17:19]}_{file[20:21]}_{file[24:27].replace("_", "")}'
        with open(f'{outdir}/{file}', 'w') as outfile:
            df.to_csv(outfile, sep = '\t')
            print(f'Pre-filtered counts of {isolate} saved to {outdir}/{file}')
        with open(f'{metadir}/{isolate}_count_distribution.tsv', 'w') as metadat:
            metadat.write('feature_type\tcount\n')
            [metadat.write(f'{k}\t{v}\n') for k, v in type_count.items()]
            print(f'RNA gene counts of {isolate} saved to {metadir}/{isolate}_count_distribution.tsv')
            
