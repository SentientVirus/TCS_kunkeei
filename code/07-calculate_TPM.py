#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:12:54 2023
A script that calculates the number of TPM based on the number of counts
@author: marina
"""

# =============================================================================
# 0. Importing packages
# =============================================================================
import os
import logging, traceback
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
# 1. Define input variables
# =============================================================================
outdir = os.path.dirname(snakemake.output[0])
outfile = snakemake.output.per_sample
mean_outfile = snakemake.output.mean
basedir = snakemake.params

if not os.path.exists(outdir):
    os.makedirs(outdir)

print(f'{outdir} exists')
# Initialize parameters for TPM calculation
rl = 76 #Average read length 
TPM = {}
    
# =============================================================================
# 2. Loop through the features in each of the files and calculate TPM
# =============================================================================
locus_tags = []
for file in sorted(snakemake.input):
    if len(locus_tags) == 0:
        check = True
    else: check = False
    with open(file) as counts:
        T = 0
        TPM_up = []
        count_df = pd.read_csv(counts, sep = '\t', skiprows = 1)
        for index, row in count_df.iterrows():
            rg = row[6] #Gene counts
            flg = row['Length']
            TPM_calc = (rg*rl)/flg
            TPM_up.append(TPM_calc*1000000)
            T += TPM_calc
            if check:
                locus_tags.append(row['Geneid'])
            if not check and row['Geneid'] not in locus_tags:
                print('OBS! Gene that is unique in this sample!')
        basename = os.path.basename(file)
        TPM_list = [TPM_val/T for TPM_val in TPM_up]
        sample = f'I{basename[17:19]}_{basename[20:21]}_{basename[24:27].replace("_", "")}'
        print(f'Calculated TPM for sample {sample}')
        TPM[sample] = TPM_list

# =============================================================================
# 3. Save TPM calculations to files
# =============================================================================
df = pd.DataFrame.from_dict(TPM, orient='index').transpose()
df['locus_tag'] = locus_tags
df = df.set_index('locus_tag')
df.to_csv(outfile, sep = '\t')
print('TPM saved to file')

# =============================================================================
# 4. Calculate mean TPMs and save to file
# =============================================================================
isol_cond = [f'{i}{c}' for i in ['01', '02', '09', '10'] for c in ['F', 'S']]
mean_df = pd.DataFrame()
for ic in isol_cond:
    underscore_ic = ic[:-1] + '_' +  ic[-1:]
    mean_df[ic] = df.filter(regex=underscore_ic).mean(axis=1)
    mean_df.to_csv(mean_outfile, sep = '\t')
print('Mean values saved to files')