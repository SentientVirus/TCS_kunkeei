#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:32:56 2023

This is a script to filter counts that removes ribosomal genes

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Importing packages
# =============================================================================

import os
import logging, traceback
from Bio import GenBank
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
# 1. Define inputs and outputs
# =============================================================================

workdir = snakemake.params.workdir #Working directory

count_dir = os.path.dirname(snakemake.input.counts[0]) #Directory with gene counts
outdir = os.path.dirname(snakemake.output.counts[0]) #Output directory
metadir = os.path.dirname(snakemake.output.summary[0]) #Directory with metadata
gbk_file = snakemake.input.gbk #Path to reference GenBank file

#Create the two output directories if they don't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
if not os.path.exists(metadir):
    os.makedirs(metadir)

# =============================================================================
# 2. Read GenBank file to get locus tags of genes annotated as ribosomal RNAs
# or tRNAs
# =============================================================================

loctags_by_type = {'rRNA': [], 'tRNA': [], 'misc_RNA': [], 'transposon': []} #Set feature types to exclude
RNA_loctags = [] #Create empty list to store locus tags to exclude

with open(gbk_file) as handle: #Open GenBank file
    for record in GenBank.parse(handle): #Loop through records in the file
        for feature in record.features: #Loop through features in each record
            if 'RNA' in feature.key or (len(feature.qualifiers) > 5 and 'transpos' in feature.qualifiers[5].value): #If the feature is an RNA-coding gene or transposon
                loctag = feature.qualifiers[0].value.replace('"', '') #Retrieve the locus tag
                check = True #Set the boolean to True
                if 'ssrA' in loctag: #If the gene is SsrA
                    loctag = feature.qualifiers[1].value.replace('"', '') #Get the locus tag
                    type_seq = 'tRNA' #Set gene type to tRNA
                    check = False #Set the check to False
                    
                print(f'{loctag} in chromosome {record.accession[0]} appended to list.', end = ' ')
                RNA_loctags.append(loctag) #Add to the list of RNA genes
                
                if check and 'rRNA' in feature.key: #If the gene is annotated as an rRNA gene
                    type_seq = 'rRNA' #Set it as the type of the feature
                elif check and 'tRNA' in feature.key: #Same for tRNA-coding genes
                    type_seq = 'tRNA'
                elif check and ('RNA' in feature.key or 'RNA' in feature.qualifiers[1].value): #Same for other RNA-coding genes
                    type_seq = 'misc_RNA'
                elif check: #If check is true, but it is not an RNA
                    type_seq = 'transposon' #Label it as a transposon
                    
                print(f'Element is a {type_seq}.')
                loctags_by_type[type_seq].append(loctag) #Add locus tag to a list of locus tags of the same type
                
        
# =============================================================================
# 3. Read count file and filter out the RNA genes, create filtered count files
# and files with the data of RNA genes 
# =============================================================================

type_count = {} #Create empty dictionary
for file in os.listdir(count_dir): #Loop through count files
    file_path = f'{count_dir}/{file}' #Retrieve the complete path to the files
    if file.endswith('featureCounts'): #If the file is part of the featureCounts output
        with open(f'{file_path}') as counts: #Open the file
            df = pd.read_csv(counts, sep = '\t', header = 1, index_col = 0) #Read it as a dataframe
            total_reads = df.iloc[:, -1:].sum().iloc[0] #Get the total number of counts
            for seqtype in ['rRNA', 'tRNA', 'misc_RNA', 'transposon']: #Loop through the types of sequences to exclude
                type_count[seqtype] = df[df.index.isin(loctags_by_type[seqtype])].iloc[:, -1:].sum().iloc[0] #Get the total number of counts
            df = df[~df.index.isin(RNA_loctags)] #Get a daframe where the genes to exclude have been filtered out
            
        isolate = f'I{file[17:19]}_{file[20:21]}_{file[24:27].replace("_", "")}' #Get the isolate number
        with open(f'{outdir}/{file}', 'w') as outfile: #Open the output file
            df.to_csv(outfile, sep = '\t') #Save the filtered dataframe to a TSV file
            print(f'Pre-filtered counts of {isolate} saved to {outdir}/{file}')
        with open(f'{metadir}/{isolate}_count_distribution.tsv', 'w') as metadat: #Open the metadata output file
            metadat.write('feature_type\tcount\n') #Write headers
            [metadat.write(f'{k}\t{v}\n') for k, v in type_count.items()] #Write the counts for each type of feature
            print(f'RNA gene counts of {isolate} saved to {metadir}/{isolate}_count_distribution.tsv')
            
