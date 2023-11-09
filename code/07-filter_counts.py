#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 10:56:25 2022

@author: marina
"""

#Script that, taking featureCounts output as input, filters low-count reads and
#builds a dataframe with the strain names as columns, the three locus tags
#as row name and the count as values.

import pandas as pd
import os
import numpy as np

#workdir = os.getcwd()
indir = snakemake.params.countdir
outdir = os.path.dirname(snakemake.output.avg_nofilter)
indir2 = snakemake.params.filtered_countdir
outdir2 = os.path.dirname(snakemake.output.avg_filtered)
ref_countfile = snakemake.input.counts_filtered[0]
#os.chdir(indir) #featureCounts_reverse/nofilter

print('Loaded input')

[os.makedirs(outd) for outd in [outdir, outdir2] if not os.path.exists(outd)]

# =============================================================================
# Here, replicate names are saved to a dictionary.
# =============================================================================
df_dict = {}
mean_fragment_length = {}
locus_tags = []
for filename in os.listdir(indir):
    if 'H3B1-04J' in filename and filename.endswith('featureCounts'):
        with open(f'{indir}/{filename}', 'r') as countfile:
            sname = filename.split('-')[4] + '_' + filename[20] + '_' + filename.split('-')[6].split('_')[1]
            features = pd.read_csv(countfile, sep = '\t', header = 1)
            locus_tags.append(features['Geneid'])
            # features = features.drop(['Start', 'End', 'Geneid'], axis = 1)
            for col in features.columns:
                if 'results' in col:
                    features.rename(columns = {col:'Results'}, inplace = True)
            df_dict[sname] = features


isolate_list = [isolate.split('_')[0] for isolate in df_dict.keys()]
isolate_list = list(np.unique(isolate_list))
conditions = ['F', 'S']

with open(snakemake.input.TPM) as tab_TPM:
    TPM_samples = pd.read_table(tab_TPM)
    print('Read TPM')
    
# =============================================================================
# First I should filter genes by TPM. We don't want TPM, but filtered counts.  
# =============================================================================

isolate_dict = {}
geneid_list = []
chr_list = []
start_list = []
end_list = []
length_list = []
results_list = []
type_dict = {}
name_dict = {}
include_dict = {}
for isolate in isolate_list:
    if isolate == '01' or isolate == '02':
        sample_type = 'mucoid'
    else:
        sample_type = 'inhibitor'
    name_dict[isolate] = sample_type
    for condition in conditions:
        samples = [f'I{n}' for n in df_dict.keys() if f'{isolate}_{condition}' in n]
        isolate_dict[f'{isolate}{condition}'] = pd.DataFrame(columns = ['Geneid', 'Chr', 'Start', 'End', 'Length'] + samples + ['Mean_results'])
        for name in df_dict.keys():
            div_part = sum(df_dict[name]['Results'])
            if f'{isolate}_{condition}' in name:
                tag = f'{isolate}{condition}'
                for index, gene in df_dict[name].iterrows():
                    x = gene['Results']
                    CPM = (x/div_part)*1000000
                    TPM = TPM_samples[TPM_samples['locus_tag'] == gene['Geneid']][tag][index]
                    if CPM < 2:
                        check = False
                        print(CPM, TPM)
                    else:
                        check = True
                    if gene['Geneid'] not in include_dict.keys():
                        include_dict[gene['Geneid']] = [check]
                    else:
                        include_dict[gene['Geneid']].append(check)

                    if gene['Geneid'] not in geneid_list:
                        geneid_list.append(gene['Geneid'])
                        chr_list.append(gene['Chr'])
                        start_list.append(gene['Start'])
                        end_list.append(gene['End'])
                        length_list.append(gene['Length'])
                        results_list.append(gene['Results'])

                isolate_dict[tag]['Geneid'] = geneid_list
                isolate_dict[tag]['Chr'] = chr_list
                isolate_dict[tag]['Start'] = start_list
                isolate_dict[tag]['End'] = end_list
                isolate_dict[tag]['Length'] = length_list
                isolate_dict[tag][f'I{name}'] = results_list
                type_dict[f'I{name}'] = isolate
                geneid_list = []
                chr_list = []
                start_list = []
                end_list = []
                length_list = []
                results_list = []
                # count_df_dict[name] = count_df
                # isolate_dict[isolate].to_csv(f'H3B1-04J_{isolate}_raw_counts.tsv', sep = '\t')
            print(name)

# =============================================================================
# Writing mean counts
# =============================================================================

gene_dict = {}
input_meta = snakemake.input.meta
new_file = snakemake.output.avg_nofilter
new_file2 = snakemake.output.avg_filtered
meta_file = snakemake.output.meta
pos_file = snakemake.output.pos
with open(ref_countfile) as keep:
    keep_row = pd.read_csv(keep, sep = '\t', index_col = 0)
    to_keep = list(keep_row.index)
    print(to_keep)
    
for i in isolate_dict.keys():
    for condition in conditions:
        isolate_dict[i]['Mean_results'] = (isolate_dict[i][isolate_dict[i].columns[5:]].sum(axis=1))/len(isolate_dict[i].columns[5:])
        for gene in isolate_dict[i]['Geneid']:
            if gene not in gene_dict.keys():
                gene_dict[gene] = list(isolate_dict[i]['Mean_results'][isolate_dict[i]['Geneid'] == gene])
            else:
                gene_dict[gene].append(float(isolate_dict[i]['Mean_results'][isolate_dict[i]['Geneid'] == gene]))
            if sum(include_dict[gene])/len(include_dict[gene]) < 0.5:
                print(gene)
                isolate_dict[i].drop(isolate_dict[i][isolate_dict[i]['Geneid'] == gene].index, inplace = True)

        counts = pd.concat([isolate_dict[i]['Geneid'], isolate_dict[i].iloc[:, 5:-1]], axis = 1)
        counts.to_csv(f'{outdir}/H3B1-04J_{i}_counts.tsv', sep = '\t', index = False)
        filtered_counts = counts[counts['Geneid'].isin(to_keep)]
        filtered_counts.to_csv(f'{outdir2}/H3B1-04J_{i}_counts.tsv', sep = '\t', index = False)
      
with open(input_meta, 'r') as exp_meta:
    df_c3 = pd.read_csv(exp_meta, sep = '\t')
    condition3 = {f'I{j}': row['batch (date, YY.MM.DD)'] for index, row in df_c3.iterrows() for j in df_dict.keys() if int(j[6:]) == row['SampleNo']}

loop_list = ['_'.join(key.split('_')[:-1]) for key in type_dict.keys()]
with open(meta_file, 'w') as metadata:
    metadata.write('sample\tcondition1\tcondition2\tcondition3\n')
    [metadata.write(f'{isol}\t{name_dict[isol[1:3]]}\t{isol[4]}\t{condition3[isol]}\n') for isol in type_dict.keys()]

#Save metadata to one file (gene name + start + end + length) and gene name +
#gene counts to a separate file    

with open(new_file, 'w') as avg_counts, open(new_file2, 'w') as avg_fc:
    avg_counts.write('\t' + '\t'.join(list(isolate_dict.keys())) + '\n')
    avg_fc.write('\t' + '\t'.join(list(isolate_dict.keys())) + '\n')
#with open(new_file, 'a') as avg_counts:
    [avg_counts.write(f'{gene}\t' + '\t'.join([str(gcount) for gcount in gene_dict[gene]]) + '\n') for gene in gene_dict.keys()]
    [avg_fc.write(f'{gene}\t' + '\t'.join([str(gcount) for gcount in gene_dict[gene]]) + '\n') for gene in gene_dict.keys() if gene in to_keep]


#Save positions to file
with open(pos_file, 'w') as pos:
    pos.write('')

with open(pos_file, 'a') as f:
    key = list(isolate_dict.keys())[0]
    df = isolate_dict[key].iloc[:, :4]
    df.set_index('Geneid', inplace = True)
    #print(df)
    df.to_csv(path_or_buf=pos_file, sep='\t')
