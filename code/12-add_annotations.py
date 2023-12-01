#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:52:12 2023

@author: marina

This script reads a GenBank file using Biopython and, taking the locus tags
in the output from DESeq2 with pandas (csv format), it adds an annotation to 
thelocus tag.

"""

# =============================================================================
# Import packages
# =============================================================================
from Bio import GenBank as gbk
import pandas as pd
import logging, traceback

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
# 1. Get inputs from Snakemake
# =============================================================================
dif_expr = snakemake.input.dif_expr
gbk_file = snakemake.input.gbk

# =============================================================================
# 2. Read the dataframes with DE results and create a list with all the
# unique locus tags
# =============================================================================
locus_tags = []
df_dict = {}
for file in dif_expr:
    print(f'Processing file {file}')
    with open(file) as csvfile:
        df = pd.read_csv(csvfile, index_col = 0)
        locus_tags = list(pd.unique(list(df.index) + locus_tags))
        df_dict[file] = df
locus_tags = sorted(locus_tags)

# =============================================================================
# 3. Create a dictionary that matches an annotation from the GenBank file
# (if it exists) to each locus tag
# =============================================================================
tag_dict = {}
with open(gbk_file) as gbk_read:
    for record in gbk.parse(gbk_read):
        for feature in record.features:
            
            if len(feature.qualifiers) > 1:
                
                # Check if the feature has a three-letter code
                annotation = ''
                
                # If it does, add it together with the protein annotation
                if 'AKUH3B104J' not in feature.qualifiers[0].value.strip('"'):
                    loctag = feature.qualifiers[1].value.strip('"')
                    gene_name = feature.qualifiers[0].value.strip('"')
                    if loctag in locus_tags:
                        for qual in feature.qualifiers:
                            if 'product' in qual.key:
                                annotation = qual.value.strip('"')
                        tag_dict[loctag] = (gene_name, annotation)
                        print(f'Found gene {gene_name}')
                        
                # Otherwise add a - symbol and the annotation
                else:
                    for qual in feature.qualifiers:
                        if 'product' in qual.key:
                            annotation = qual.value.strip('"')
                    loctag = feature.qualifiers[0].value.strip('"')
                    tag_dict[loctag] = ('-', annotation)
                    print(f'Added hypothetical protein {loctag}')
                        
# =============================================================================
# 4. Add annotations to the locus tags in the dataframes and 
# save to tsv files
# =============================================================================

for key in df_dict.keys():
    df = df_dict[key]
    new_column = []
    locus_tags = list(df.index)
    new_list = [tag_dict[loctag][0] for loctag in locus_tags]
    annot = [tag_dict[loctag][1] for loctag in locus_tags]
    df.insert(len(df.columns), 'gene_name', new_list, True)
    df.insert(len(df.columns), 'annotation', annot, True)
    outfile = key.replace('.csv', '_annotated.tsv')
    df.to_csv(outfile, sep='\t')
    print(f'Saved annotations for comparison {key.replace(".csv", "")}')
