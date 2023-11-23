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
import os

# =============================================================================
# Section to define inputs
# =============================================================================
home = os.path.expanduser('~')
workdir = home + '/snpseq00064'
csv_dir = workdir + '/results/DE'
gbk_file = home + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff' 

# =============================================================================
# Section to include logging
# =============================================================================

# =============================================================================
# Function to update dataframes with annotations
# =============================================================================

# =============================================================================
# 1. Extract the dataframes with DE results and create a list with all the
# unique locus tags
# =============================================================================
locus_tags = []
df_dict = {}
for file in os.listdir(csv_dir):
    if file.endswith('.csv'):
        print(file)
        with open(f'{csv_dir}/{file}') as csvfile:
            df = pd.read_csv(csvfile, index_col = 0)
            locus_tags = list(pd.unique(list(df.index) + locus_tags))
            df_dict[file] = df
locus_tags = sorted(locus_tags)

# =============================================================================
# 2. Create a dictionary that matches an annotation from the GenBank file
# (if it exists) to each locus tag
# =============================================================================
tag_dict = {}
with open(gbk_file) as gbk_read:
    for record in gbk.parse(gbk_read):
        for feature in record.features:
            if len(feature.qualifiers) > 1:
                if 'AKUH3B104J' not in feature.qualifiers[0].value.strip('"'):
                    loctag = feature.qualifiers[1].value.strip('"')
                    if loctag in locus_tags:
                        tag_dict[loctag] = feature.qualifiers[0].value.strip('"')
                else:
                    loctag = feature.qualifiers[0].value.strip('"')
                    tag_dict[loctag] = ''
                        
# =============================================================================
# 3. Use function to add annotations to the locus tags in the dataframes
# =============================================================================
for key in df_dict.keys():
    df = df_dict[key]
    new_column = []
    locus_tags = list(df.index)
    new_list = [tag_dict[loctag] for loctag in locus_tags]
    df.insert(len(df.columns), 'gene_name', new_list, True)
    outfile = key.replace('.csv', '_annotated.tsv')
    df.to_csv(f'{csv_dir}/{outfile}', sep='\t') 

# =============================================================================
# 4. Save to .tsv files
# =============================================================================

# =============================================================================
# Do the same with count files (create files with locus tags, counts and 
# annotation)
# =============================================================================
