#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:12:54 2023
A script that calculates the number of TPM based on the number of counts
@author: marina
"""

# =============================================================================
# A Python script to calculate TPM
# =============================================================================
import os
import pandas as pd

# =============================================================================
# 1. Define input variables
# =============================================================================
#workdir = '../featureCounts_reverse/nofilter' #Change to snakemake input directory
outdir = os.path.dirname(snakemake.output[0]) #'../results/TPM' #Change to snakemake output directory
outfile = snakemake.output.per_sample #'TPM_per_sample.tsv'
mean_outfile = snakemake.output.mean #'mean_TPM.tsv'
basedir = snakemake.params

if not os.path.exists(outdir):
    os.makedirs(outdir)

# Initialize parameters for TPM calculation
rl = 76 #Average read length 
TPM = {}
#sum_dict = {}

# # =============================================================================
# # 2. Open the summary files to get the total number of reads per sample
# # =============================================================================

# os.chdir(basedir)
# for file in sorted(snakemake.input.summary):
#     #if file.endswith('summary'): #This line should be removed when I implement snakemake
#     with open(file) as summary:
#         sample = file.replace('.summary', '') #f'I{file[17:19]}_{file[20:21]}_{file[24:27].replace("_", "")}'
#         sum_stats = pd.read_csv(summary, sep = '\t').iloc[:,1].sum()
#         sum_dict[sample] = sum_stats
#         # if sample == 'I01_F_S5':
#         #     print(sum_stats, sum_stats.sum())
#         print(f'Read {sample}')
    
# =============================================================================
# 2. Loop through the features in each of the files and calculate TPM
# =============================================================================
locus_tags = []
for file in sorted(snakemake.input):
    #if file.endswith('featureCounts'): #This line should be removed when I implement snakemake
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

# =============================================================================
# 4. Calculate mean TPMs and save to file
# =============================================================================
isol_cond = [f'{i}{c}' for i in ['01', '02', '09', '10'] for c in ['F', 'S']]
mean_df = pd.DataFrame()
for ic in isol_cond:
    underscore_ic = ic[:-1] + '_' +  ic[-1:]
    mean_df[ic] = df.filter(regex=underscore_ic).mean(axis=1)
    mean_df.to_csv(mean_outfile, sep = '\t')
print('Data saved to files')