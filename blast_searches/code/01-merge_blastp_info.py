#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 16:14:23 2025

This script takes two files from BLAST, one with alignment hits and another
from hit descriptions, and it merges the information together.
Environment: alignment_tree.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, pandas as pd

# =============================================================================
# 1. Define input files
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/blast_searches'
hitfile = f'{workdir}/blast_RR-TF.csv'
#hitfile0 = f'{workdir}/GS1_H1B3-02M_Lactobacillaceae.csv'
#extra_hitfiles = [hitfile0] 
descfile = f'{workdir}/blast_RR-TF_descriptions.csv'
#descfile0 = f'{workdir}/GS1_H1B3-02M_hit_descriptions.csv'
#extra_descfiles = [descfile0]
outfile = f'{workdir}/blastp_RR-TF_results.csv'

# =============================================================================
# 2. Read the files and merge the information in both of them
# =============================================================================

with open(hitfile) as hits:
    hit_df = pd.read_csv(hitfile, header = None, sep = ',')
    hit_df.rename(columns={0: 'Query ID', 1: 'Subject accession', 
                           2: '% Identity', 3: 'Align. len', 4: 'Mismatches', 
                           5: 'No. gaps', 6: 'Query start', 7: 'Query end', 
                           8: 'Subject start', 9: 'Subject end', 10: 'E-value',
                           11: 'Bit score', 12: '% Positives'
                           }, inplace=True)
    hit_df.drop(columns='Query ID', inplace=True)
    
with open(descfile) as desc:
    desc_df = pd.read_csv(descfile, sep = ',')
    desc_df.rename(columns={'Accession  ': 'Subject accession', 
                            'Query Cover': 'Query cover',
                            'Max Score': 'Max. score',
                            'Total Score': 'Total score',
                            'Acc. Len': 'Acc. len',
                            'Scientific Name': 'Scientific name'}, inplace=True)
    desc_df.drop(columns=['Description', 'Per. ident', 'E value'], inplace=True)
    desc_df['Subject accession'] = desc_df['Subject accession'].apply(lambda x: x.split('/')[4].split('?')[0])
    
df = hit_df.merge(desc_df, on='Subject accession', how='left')
    
# extra_hit_dfs = []
# for i in range(len(extra_hitfiles)):
#     file = extra_hitfiles[i]
#     desc_file = extra_descfiles[i]
    
#     extra_df = pd.read_csv(file, header = None, sep = ',')
#     extra_df.rename(columns={0: 'Query ID', 1: 'Subject accession',
#                            2: '% Identity', 3: 'Align. len', 4: 'Mismatches',
#                            5: 'No. gaps', 6: 'Query start', 7: 'Query end',
#                            8: 'Subject start', 9: 'Subject end', 10: 'E-value',
#                            11: 'Bit score', 12: '% Positives'
#                            }, inplace=True)
    
#     add_df = extra_df[~extra_df['Subject accession'].isin(df['Subject accession'])]
    
#     if len(add_df) > 0:
#         with open(desc_file) as desc:
#             desc_extra = pd.read_csv(desc, sep = ',')
#             desc_extra.rename(columns={'Accession  ': 'Subject accession', 
#                                     'Query Cover': 'Query cover',
#                                     'Max Score': 'Max. score',
#                                     'Total Score': 'Total score',
#                                     'Acc. Len': 'Acc. len',
#                                     'Scientific Name': 'Scientific name'}, inplace=True)
#             desc_extra.drop(columns=['Description', 'Per. ident', 'E value'], inplace=True)
#             desc_extra['Subject accession'] =  desc_extra['Subject accession'].apply(lambda x: x.split('/')[4].split('?')[0])
#         new_df = add_df.merge(desc_extra, on='Subject accession', how='left')
#         df = pd.concat([df, new_df], axis=0)
    


reorder = ['Subject accession', 'Scientific name', '% Identity', 
           'Query cover', 'Acc. len', 'Align. len', 'Mismatches', 'No. gaps', 
           'Query start', 'Query end', 'Subject start', 'Subject end', 
           'E-value', 'Bit score', 'Max. score', 'Total score', '% Positives']

df = df.reindex(columns = reorder)

with open(outfile, 'w') as out_csv:
    df.to_csv(out_csv, sep = ',', index = False)
