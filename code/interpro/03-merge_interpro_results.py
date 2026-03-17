#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 17:03:59 2026

Script to merge the InterProScan results for each file into a single file.
The script is numbered as "03" to remind the user that there is an intermediate
step, which is running InterProScan.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
import pandas as pd

# =============================================================================
# 2. Define paths to inputs and outputs
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/snpseq00064'

indir = f'{workdir}/results/InterProScan'
infiles = [f'{indir}/{file}' for file in os.listdir(indir) if 'proteins' in file and file.endswith('.tsv') and 'complete' not in file]
infiles = sorted(infiles)
outfile = f'{indir}/H3B1-04J_complete_proteins.tsv'

# =============================================================================
# 3. Create a list with the desired column names
# =============================================================================

col_list = ['Protein accession', 'Sequence MD5 digest', 'Sequence length',
              'Analysis', 'Signature accession', 'Signature description',
              'Start location', 'Stop location', 'Score', 'Status', 'Date',
              'InterPro accession', 'InterPro description', 'GO annotations',
              'Pathways annotations']

# =============================================================================
# 4. Merge InterProScan results into a single dataframe
# =============================================================================

df_list = [pd.read_csv(file, sep = '\t', header = None) for file in infiles] #Read the files as dataframes and save them to a list

df = pd.concat(df_list, ignore_index = True, sort = True, axis = 0) #Merge the dataframes
df.columns = col_list #Change the column names

#Split the protein ID column into two, one with the protein ID and another with the locus tag
df[['Protein accession', 'Locus tag']] = df['Protein accession'].str.split('/', n = 1, expand = True)

df_final = df[[col_list[0], 'Locus tag'] + col_list[1:]] #Move the locus tag column to the position after the protein ID column

df_final.sort_values(df_final.columns[1], axis = 0, inplace = True, ignore_index = True) #Sort the dataframe by locus tag

# =============================================================================
# 5. Write the dataframe to an output file
# =============================================================================

with open(outfile, 'w') as handle: #Open the output file in write mode
    df_final.to_csv(handle, sep ='\t', index = False) #Write the dataframe to the file without the index column