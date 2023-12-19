#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 10:56:25 2022

Script that, taking featureCounts output as input, filters low-count reads and
builds a dataframe with the strain names as columns, the three locus tags
as row name and the count as values.

@author: marina
"""

# =============================================================================
# 0. Importing packages
# =============================================================================

import pandas as pd
import os
import numpy as np
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

indir = snakemake.params.countdir
outdir = os.path.dirname(snakemake.output.avg_nofilter)
indir2 = snakemake.params.filtered_countdir
outdir2 = os.path.dirname(snakemake.output.avg_filtered)
ref_countfile = snakemake.input.counts_filtered[0]

print('Loaded input')

[os.makedirs(outd) for outd in [outdir, outdir2] if not os.path.exists(outd)]

# =============================================================================
# 2. Save data to a dictionary under replicate names
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
            for col in features.columns:
                if 'results' in col:
                    features.rename(columns = {col:'Results'}, inplace = True)
            df_dict[sname] = features
            print(f'Processed {indir}/{filename}')


isolate_list = [isolate.split('_')[0] for isolate in df_dict.keys()]
isolate_list = list(np.unique(isolate_list))
conditions = ['F', 'S']

# =============================================================================
# 3. Read TPM data
# =============================================================================
with open(snakemake.input.TPM) as tab_TPM:
    TPM_samples = pd.read_table(tab_TPM)
    print('Read TPM')
    
# =============================================================================
# 4. Filter counts by TPM  
# =============================================================================

# Variable definition
isolate_dict = {}
# geneid_list = []
# chr_list = []
# start_list = []
# end_list = []
# length_list = []
# results_list = []
type_dict = {}
name_dict = {}
include_dict = {}

# Loop through isolates and conditions
for isolate in isolate_list:
    if isolate == '01' or isolate == '02':
        sample_type = 'mucoid'
    else:
        sample_type = 'inhibitor'
    name_dict[isolate] = sample_type
    for condition in conditions:
        samples = [f'I{n}' for n in df_dict.keys() if f'{isolate}_{condition}' in n]
        isolate_dict[f'{isolate}{condition}'] = pd.DataFrame(columns = ['Geneid', 'Chr', 'Start', 'End', 'Length'] + samples + ['Mean_results'])
        
        # Filter counts by TPM
        for name in df_dict.keys():
            #div_part = sum(df_dict[name]['Results'])
            
            # Create variables to be used
            if f'{isolate}_{condition}' in name:
                tag = f'{isolate}{condition}'
                geneid_list = []
                chr_list = []
                start_list = []
                end_list = []
                length_list = []
                results_list = []
                for index, gene in df_dict[name].iterrows():
                    #x = gene['Results']
                    #CPM = (x/div_part)*1000000
                    
                    # For each sample, check if the TPM criterion is met
                    TPM = TPM_samples[TPM_samples['locus_tag'] == gene['Geneid']][tag][index]
                    if TPM < 10:
                        check = False
                    else:
                        check = True
                        
                    # Create list to check how many isolates fulfil the condition
                    if gene['Geneid'] not in include_dict.keys():
                        include_dict[gene['Geneid']] = [check]
                    else:
                        include_dict[gene['Geneid']].append(check)
                    
                    # Add gene to the list of genes to consider if it is not there already
                    if gene['Geneid'] not in geneid_list:
                        geneid_list.append(gene['Geneid'])
                        chr_list.append(gene['Chr'])
                        start_list.append(gene['Start'])
                        end_list.append(gene['End'])
                        length_list.append(gene['Length'])
                        results_list.append(gene['Results'])
                
                # Save gene information to a dictionary
                isolate_dict[tag]['Geneid'] = geneid_list
                isolate_dict[tag]['Chr'] = chr_list
                isolate_dict[tag]['Start'] = start_list
                isolate_dict[tag]['End'] = end_list
                isolate_dict[tag]['Length'] = length_list
                isolate_dict[tag][f'I{name}'] = results_list
                type_dict[f'I{name}'] = isolate
                
                # # Empty variables
                # geneid_list = []
                # chr_list = []
                # start_list = []
                # end_list = []
                # length_list = []
                # results_list = []
            print(f'Filtered counts of {name}')

# =============================================================================
# 5. Write counts per isolate
# =============================================================================

# Variable definitions
gene_dict = {}
input_meta = snakemake.input.meta
new_file = snakemake.output.avg_nofilter
new_file2 = snakemake.output.avg_filtered
meta_file = snakemake.output.meta
pos_file = snakemake.output.pos

# Read from files which genes should be kept
with open(ref_countfile) as keep:
    keep_row = pd.read_csv(keep, sep = '\t', index_col = 0)
    to_keep = list(keep_row.index)
    
# Loop through isolates and conditions to calculate mean counts    
for i in isolate_dict.keys():
    isolate_dict[i]['Mean_results'] = (isolate_dict[i][isolate_dict[i].columns[5:]].sum(axis=1))/len(isolate_dict[i].columns[5:])

    # Loop through genes to retrieve mean counts for that isolate
    for gene in isolate_dict[i]['Geneid']:
        if gene not in gene_dict.keys():
            gene_dict[gene] = [float(isolate_dict[i]['Mean_results'][isolate_dict[i]['Geneid'] == gene])]
        else:
            gene_dict[gene].append(float(isolate_dict[i]['Mean_results'][isolate_dict[i]['Geneid'] == gene]))
            
        # Remove genes if they don't meet the TPM threshold in >50% samples
        if sum(include_dict[gene])/len(include_dict[gene]) < 0.5:
            print(f'Excluded {gene} in {i} due to low TPM')
            isolate_dict[i].drop(isolate_dict[i][isolate_dict[i]['Geneid'] == gene].index, inplace = True)
    
    # For each isolate and condition, save the counts of the five replicates to a file
    counts = pd.concat([isolate_dict[i]['Geneid'], isolate_dict[i].iloc[:, 5:-1]], axis = 1)
    counts.to_csv(f'{outdir}/H3B1-04J_{i}_counts.tsv', sep = '\t', index = False)
    filtered_counts = counts[counts['Geneid'].isin(to_keep)]
    filtered_counts.to_csv(f'{outdir2}/H3B1-04J_{i}_counts.tsv', sep = '\t', index = False)
    print(f'Saved counts to H3B1-04J_{i}_counts.tsv')
      
# =============================================================================
# 6. Create metadata file
# =============================================================================

# Retrieve infomation from file
with open(input_meta, 'r') as exp_meta:
    df_c3 = pd.read_csv(exp_meta, sep = '\t')
    condition3 = {f'I{j}': row['batch (date, YY.MM.DD)'] for index, row in df_c3.iterrows() for j in df_dict.keys() if int(j[6:]) == row['SampleNo']}

# Save information to a new file
loop_list = ['_'.join(key.split('_')[:-1]) for key in type_dict.keys()]
with open(meta_file, 'w') as metadata:
    metadata.write('sample\tcondition1\tcondition2\tcondition3\n')
    [metadata.write(f'{isol}\t{name_dict[isol[1:3]]}\t{isol[4]}\t{condition3[isol]}\n') for isol in type_dict.keys()]
    print('Created metadata file')

# =============================================================================
# 7. Save mean counts (collapsing replicates) into a single file   
# =============================================================================

with open(new_file, 'w') as avg_counts, open(new_file2, 'w') as avg_fc:
    print(isolate_dict)
    avg_counts.write('\t' + '\t'.join(list(isolate_dict.keys())) + '\n')
    avg_fc.write('\t' + '\t'.join(list(isolate_dict.keys())) + '\n')
    [avg_counts.write(f'{gene}\t' + '\t'.join([str(gcount) for gcount in gene_dict[gene]]) + '\n') for gene in gene_dict.keys()]
    [avg_fc.write(f'{gene}\t' + '\t'.join([str(gcount) for gcount in gene_dict[gene]]) + '\n') for gene in gene_dict.keys() if gene in to_keep]
    print('Created files with average gene counts per isolate and condition')


# =============================================================================
# 8. Save extra metadata to be used as input for the saturation analysis
#                   (gene name + start + end + length)
# =============================================================================

# Create empty file
with open(pos_file, 'w') as pos:
    pos.write('')

# Write information to file
with open(pos_file, 'a') as f:
    key = list(isolate_dict.keys())[0]
    df = isolate_dict[key].iloc[:, :4]
    df.set_index('Geneid', inplace = True)
    df.to_csv(path_or_buf=pos_file, sep='\t')
    print('Created detailed metadata file for the saturation analysis')