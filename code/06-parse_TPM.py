#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 15:20:38 2023

@author: marina
"""
import os
import pandas as pd
from statistics import mean

indir = 'results/TPM'
outfile = 'mean_TPM.tsv'

samples = ['01', '02', '09', '10']
conditions = ['F', 'S']
sc = [sample + condition for sample in samples for condition in conditions]
tpms = {}

for smpl in sc:
    tpms[smpl] = {}

for file in os.listdir(indir):
    if file.endswith('.out'):
        with open(f'{indir}/{file}') as TPMs:
            s_type = file[17:19] + file[20]
            df = pd.read_csv(TPMs, sep = '\t')
            df.set_index(df.columns[0], inplace=True)
            df.index.name = None
            if tpms[s_type] == {}:
                tpms[s_type] = df['UniqueTPM'].to_dict()
                tpms[s_type] = {k: [v] for k, v in tpms[s_type].items()}
            else: [tpms[s_type][index].append(df['UniqueTPM'][index]) for index in list(df.index)]

# with open(f'{indir}/{outfile}', 'w') as TPM_summary:
#     TPM_summary.write('\t'.join(sc) + '\n') 
            
dict_list = []
sample_dict = {}
i = 0
# Now I only have to calculate the mean TPMs per type and save this to a tab file
for stype in tpms.keys():
    new_dict = {}
    for key in tpms[stype].keys():
        new_dict[key] = mean(tpms[stype][key])
    dict_list.append(new_dict)
        
TPM_df = pd.DataFrame(dict_list).T
TPM_df.columns = sc
TPM_df.index.name = 'locus_tag'

TPM_df.to_csv(f'{indir}/{outfile}', sep = '\t', index = True, lineterminator = '\n')
