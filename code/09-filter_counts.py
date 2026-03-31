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

indir = snakemake.params.countdir #Input directory
outdir = os.path.dirname(snakemake.output.avg_nofilter) #Output directory
indir2 = snakemake.params.filtered_countdir #Second input directory
outdir2 = os.path.dirname(snakemake.output.avg_filtered) #Second output directory
ref_countfile = snakemake.input.counts_filtered[0] #Reference file with gene counts

print('Loaded input')

[os.makedirs(outd) for outd in [outdir, outdir2] if not os.path.exists(outd)] #Create output directories if they don't exist

TPM_threshold = 10 #Threshold below which genes are filtered out if TPM < threshold in >50% samples

# =============================================================================
# 2. Save data to a dictionary under replicate names
# =============================================================================

df_dict = {} #Create an empty dictionary to store dataframes
locus_tags = [] #Create an empty list to store locus tags

for filename in os.listdir(indir): #Loop through input files
    if 'H3B1-04J' in filename and filename.endswith('featureCounts'): #If the file contains count data
        with open(f'{indir}/{filename}', 'r') as countfile: #Open the file
            sname = filename.split('-')[4] + '_' + filename[20] + '_' + filename.split('-')[6].split('_')[1] #Retrieve the sample name
            features = pd.read_csv(countfile, sep = '\t', header = 1) #Read the file as a dataframe
            locus_tags.append(features['Geneid']) #Add the locus tags of the genes to a list
            for col in features.columns: #Loop through the columns in the dataframe
                if 'results' in col: #If it is the "results" column, rename it
                    features.rename(columns = {col:'Results'}, inplace = True)
            df_dict[sname] = features #Save the dataframe to the dictionary
            print(f'Processed {indir}/{filename}')


isolate_list = [isolate.split('_')[0] for isolate in df_dict.keys()] #Get a list of isolates
isolate_list = list(np.unique(isolate_list)) #Keep only one occurrence of each
conditions = ['F', 'S'] #Create a list of conditions

# =============================================================================
# 3. Read TPM data
# =============================================================================

with open(snakemake.input.TPM) as tab_TPM: #Open TPM file
    TPM_samples = pd.read_table(tab_TPM) #Read as dataframe
    print('Read TPM')
    
# =============================================================================
# 4. Filter counts by TPM  
# =============================================================================

# Create empty dictionaries to store information
isolate_dict = {}
type_dict = {}
name_dict = {}
include_dict = {}

# Loop through isolates and conditions
for isolate in isolate_list:
    if isolate == '01' or isolate == '02': #Depending on isolate name, assign phenotype
        sample_type = 'mucoid'
    else:
        sample_type = 'inhibitor'
    name_dict[isolate] = sample_type #Store phenotype information in a dictionary
    
    for condition in conditions: #Loop through conditions
        samples = [f'I{n}' for n in df_dict.keys() if f'{isolate}_{condition}' in n] #Retrieve sample names
        isolate_dict[f'{isolate}{condition}'] = pd.DataFrame(columns = ['Geneid', 'Chr', 'Start', 'End', 'Length'] + samples + ['Mean_results']) #Create an empty dictionary and add sample information
        
        for name in df_dict.keys(): #Loop through keys in the first dictionary
            
            if f'{isolate}_{condition}' in name: #If the isolate and condition are in the name
                tag = f'{isolate}{condition}' #Create tag with this information
                geneid_list = [] #Create lists that will fill the columns in the isolate dictionary
                chr_list = []
                start_list = []
                end_list = []
                length_list = []
                results_list = []
                
                for index, gene in df_dict[name].iterrows(): #Loop through genes in the first dictionary
                    TPM = TPM_samples[TPM_samples['locus_tag'] == gene['Geneid']][tag][index] #Get the TPM for the gene
                    if TPM < 10: #If the TPM count is lower than 10 in that sample
                        check = False #Set boolean to False
                    else:
                        check = True #Else, set it to True
                        
                    # Create list of booleans to check how many isolates fulfill the condition
                    if gene['Geneid'] not in include_dict.keys():
                        include_dict[gene['Geneid']] = [check]
                    else:
                        include_dict[gene['Geneid']].append(check)
                    
                    # Add gene to the lists of genes to add to the dataframe if it is not there already
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

            print(f'Filtered counts of {name}')

# =============================================================================
# 5. Write counts per isolate
# =============================================================================

# Output definitions
gene_dict = {}
input_meta = snakemake.input.meta
new_file = snakemake.output.avg_nofilter
new_file2 = snakemake.output.avg_filtered
meta_file = snakemake.output.meta
pos_file = snakemake.output.pos

# Read from files which genes should be kept
with open(ref_countfile) as keep: #Open the file with pre-filtered counts
    keep_row = pd.read_csv(keep, sep = '\t', index_col = 0) #Read as dataframe
    to_keep = list(keep_row.index) #List of locus tags (index column)
        
for i in isolate_dict.keys(): # Loop through isolates and conditions
    print(i)
    print(isolate_dict[i])
    print(isolate_dict[i].columns[5:])
    isolate_dict[i]['Mean_results'] = (isolate_dict[i][isolate_dict[i].columns[5:]].sum(axis=1))/(len(isolate_dict[i].columns[5:])-1) #Calculate mean counts

    for gene in isolate_dict[i]['Geneid']: #Loop through genes
        if gene not in gene_dict.keys(): #If the gene is not in the gene dictionary
            gene_dict[gene] = [float(isolate_dict[i]['Mean_results'][isolate_dict[i]['Geneid'] == gene].iloc[0])] #Assign the mean counts to the dictionary
        else:
            gene_dict[gene].append(float(isolate_dict[i]['Mean_results'][isolate_dict[i]['Geneid'] == gene].iloc[0]))  #Else, append them to a list
            
        if sum(include_dict[gene])/len(include_dict[gene]) < 0.5: #If a gene doesn't meet the TPM threshold in >50% of the samples
            print(f'Excluded {gene} in {i} due to low TPM')
            isolate_dict[i].drop(isolate_dict[i][isolate_dict[i]['Geneid'] == gene].index, inplace = True) #Exclude it
    
    print(isolate_dict[i])
    
    # For each isolate and condition, save the counts of the five replicates to a file
    counts = pd.concat([isolate_dict[i]['Geneid'], isolate_dict[i].iloc[:, 5:-1]], axis = 1) #Create a dataframe concatenating the samples of each isolate and condition
    counts.to_csv(f'{outdir}/H3B1-04J_{i}_counts.tsv', sep = '\t', index = False) #Save the dataframe to TSV file
    
    # Same, but with the pre-filtered counts
    filtered_counts = counts[counts['Geneid'].isin(to_keep)]
    filtered_counts.to_csv(f'{outdir2}/H3B1-04J_{i}_counts.tsv', sep = '\t', index = False)
    
    print(f'Saved counts to H3B1-04J_{i}_counts.tsv')
      
# =============================================================================
# 6. Create metadata file
# =============================================================================

# Retrieve  batch date infomation from file
with open(input_meta, 'r') as exp_meta: #Open file with metadata
    df_c3 = pd.read_csv(exp_meta, sep = '\t') #Read file as a dataframe
    condition3 = {f'I{j}': row['batch (date, YY.MM.DD)'] for index, row in df_c3.iterrows() for j in df_dict.keys() if int(j[6:]) == row['SampleNo']} #Retrieve batch date and assign it to each sample

# Save information to a new file
with open(meta_file, 'w') as metadata: #Open the output metadata file
    metadata.write('sample\tcondition1\tcondition2\tcondition3\n') #Write headers
    [metadata.write(f'{isol}\t{name_dict[isol[1:3]]}\t{isol[4]}\t{condition3[isol]}\n') for isol in type_dict.keys()] #Write info of the three conditions (isolate, carbon source and date)
    print('Created metadata file')

# =============================================================================
# 7. Save mean counts (collapsing replicates) into a single file   
# =============================================================================

with open(new_file, 'w') as avg_counts, open(new_file2, 'w') as avg_fc: #Open output count files (for non-filtered and filtered counts)
    print(isolate_dict)
    # Write headers
    avg_counts.write('\t' + '\t'.join(list(isolate_dict.keys())) + '\n')
    avg_fc.write('\t' + '\t'.join(list(isolate_dict.keys())) + '\n')
    #Write mean results to files
    [avg_counts.write(f'{gene}\t' + '\t'.join([str(gcount) for gcount in gene_dict[gene]]) + '\n') for gene in gene_dict.keys()]
    [avg_fc.write(f'{gene}\t' + '\t'.join([str(gcount) for gcount in gene_dict[gene]]) + '\n') for gene in gene_dict.keys() if gene in to_keep]
    print('Created files with average gene counts per isolate and condition')


# =============================================================================
# 8. Save extra metadata to be used as input for the saturation analysis
#                   (gene name + start + end + length)
# =============================================================================

# Write information to file
with open(pos_file, 'w') as f: #Open output file in write mode
    key = list(isolate_dict.keys())[0] #Retrieve dictionary keys
    df = isolate_dict[key].iloc[:, :4] #Retrieve the columns with gene positions
    df.set_index('Geneid', inplace = True) #Use the gene ID as index
    df.to_csv(path_or_buf=pos_file, sep='\t') #Write dataframe to file
    print('Created detailed metadata file for the saturation analysis')
