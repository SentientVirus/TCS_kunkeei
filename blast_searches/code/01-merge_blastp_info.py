#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 16:14:23 2025

This script takes two files from BLAST, one with alignment hits and another
from hit descriptions, and it merges the information together.
Environment: pixi_phylo/default.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, pandas as pd
import logging, traceback
import sys

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0]

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    '''Function to handle exceptions'''
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger = logging.getLogger() #Create a logger

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ])) #Format to show exceptions

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Define input files
# =============================================================================

hitfile = snakemake.input.hits #Path to the file with BLASTp hits (query cover, % id, etc.)
descfile = snakemake.input.descriptions #Path to the file with hit descriptions (species, annotations, etc.)
outfile = snakemake.output[0] #Path to the output file

# =============================================================================
# 2. Read the files and merge the information in both of them
# =============================================================================

print(f'Reading {hitfile} as a dataframe... (1/4)')
with open(hitfile) as hits: #Open the file with hit information
    hit_df = pd.read_csv(hitfile, header = None, sep = ',') #Read the file as a dataframe
    hit_df.rename(columns={0: 'Query ID', 1: 'Subject accession', 
                           2: '% Identity', 3: 'Align. len', 4: 'Mismatches', 
                           5: 'No. gaps', 6: 'Query start', 7: 'Query end', 
                           8: 'Subject start', 9: 'Subject end', 10: 'E-value',
                           11: 'Bit score', 12: '% Positives' #Re-name the columns
                           }, inplace = True)
    #Remove the column with the query (same for all hits)
    hit_df.drop(columns ='Query ID', inplace = True)
print('Done! (1/4)')    
    
print(f'Reading {descfile} as a dataframe... (2/4)')
with open(descfile) as desc: #Open the file with the hit descriptions
    desc_df = pd.read_csv(descfile, sep = ',') #Read the file as a daraframe
    desc_df.rename(columns={'Accession  ': 'Subject accession', 
                            'Query Cover': 'Query cover',
                            'Max Score': 'Max. score',
                            'Total Score': 'Total score',
                            'Acc. Len': 'Acc. len',
                            'Scientific Name': 'Scientific name'},
                   inplace = True) #Re-name the columns
    desc_df.drop(columns = ['Description', 'Per. ident', 'E value'], 
                 inplace = True) #Remove columns not wanted in the output
    #Retrieve the accession of the hit from one of the columns
    desc_df['Subject accession'] = desc_df['Subject accession'].apply(lambda x: x.split('/')[4].split('?')[0])
print('Done! (2/4)')
    
print('Merging the two dataframes... (3/4)')
df = hit_df.merge(desc_df, on = 'Subject accession', how = 'left') #Merge the two dataframes

#Re-order dataframe columns
reorder = ['Subject accession', 'Scientific name', '% Identity', 
           'Query cover', 'Acc. len', 'Align. len', 'Mismatches', 'No. gaps', 
           'Query start', 'Query end', 'Subject start', 'Subject end', 
           'E-value', 'Bit score', 'Max. score', 'Total score', '% Positives']

df = df.reindex(columns = reorder) #Reset the index of the dataframe
print('Done! (3/4)')

print('Writing the merged dataframe to file... (4/4)')
with open(outfile, 'w') as out_csv: #Open the output file
    df.to_csv(out_csv, sep = ',', index = False) #Write the dataframe to the file without the index column
print('All done!')