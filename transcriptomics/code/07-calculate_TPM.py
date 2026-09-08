#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:12:54 2023

A script that calculates the number of TPM based on the number of counts

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Importing packages
# =============================================================================

import os, sys
import logging, traceback
import pandas as pd
from Bio import SeqIO

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

    # Create a logger
    logger = logging.getLogger()

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# 1. Define functions to be used in the script
# =============================================================================

def parse_gbk(main_gbk: str, ref_gbk: str):
    """
    Function that takes the path to a GenBank file and parses its contents to
    create a dictionary that maps gene names to locus tags.

    Parameters
    ----------
    ref_gbk : str
        Path to the input GenBank file.

    Returns
    -------
    ref_dict : dict
        A dictionary where locus tags are the keys and gene name annotations
        are the values.

    """
    
    
    ref_dict = {} #Create output dictionary
    with open(main_gbk) as gbk: #Open input GenBank
        for contig in SeqIO.parse(gbk, 'genbank'): #Loop through records (contigs) in the GenBank
            for CDS in contig.features: #Loop through features
                if (CDS.type == 'CDS' or 'RNA' in CDS.type) and 'locus_tag' in CDS.qualifiers.keys(): #If the feature has a locus tag
                    loctag = CDS.qualifiers['locus_tag'][0] #Retrieve it
                    if 'gene' in CDS.qualifiers.keys(): #If the feature has the "gene" property (gene name)
                        gname = CDS.qualifiers['gene'][0] #Retrieve it
                    else: gname = '-' #Else, set "-" as gene name
                    if 'product' in CDS.qualifiers.keys(): #If the feature has a "description" property
                        description = CDS.qualifiers['product'][0] #Retrieve it
                    else: description = '-' #Else, set "-" as gene name
                    
                    ref_dict[loctag] = [gname, description] #Assign the gene name to the locus tag in the output dictionary
    
    
    with open(ref_gbk) as gbk: #Open input GenBank
        for contig in SeqIO.parse(gbk, 'genbank'): #Loop through records (contigs) in the GenBank
            for CDS in contig.features: #Loop through features
                if CDS.type == 'CDS' and 'locus_tag' in CDS.qualifiers.keys(): #If the feature has a locus tag
                    loctag = CDS.qualifiers['locus_tag'][0] #Retrieve it
                    if 'gene' in CDS.qualifiers.keys(): #If the feature has the "gene" property (gene name)
                        ref_dict[loctag][0] = CDS.qualifiers['gene'][0] #Update the dictionary
                    if 'product' in CDS.qualifiers.keys(): #If the feature has a "description" property
                        ref_dict[loctag][1] = CDS.qualifiers['product'][0] #Update the dictionary
                    
    return ref_dict

# =============================================================================
# 2. Define input variables
# =============================================================================

outdir = os.path.dirname(snakemake.output[0]) #Output directory
outfile = snakemake.output.per_sample #Output file with per sample statistics
mean_outfile = snakemake.output.mean #Output file with per condition mean statistics
annot_outfile = snakemake.output.annot #Same as above, but with gene annotations
basedir = snakemake.params #Working directory
base_gbk = snakemake.input.base_gbk #Original GenBank file
ref_gbk = snakemake.input.gbk #Path to the GenBank file to add annotations

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

TPM = {} # Initialize dictionary for TPM calculation
    
# =============================================================================
# 3. Loop through the features in each of the files and calculate TPM
# =============================================================================
locus_tags = [] #List to store locus tags
for file in sorted(snakemake.input.counts): #Loop through input files
    if len(locus_tags) == 0: #If the list with locus tags is empty
        check = True #Set the boolean to True
    else: check = False #Else, set it to False
    
    with open(file) as counts: #Open input file
        count_df = pd.read_csv(counts, sep = '\t', skiprows = 1) #Read the file as a dataframe
        T = 0 #Variable with the sum of normalized counts
        TPM_up = [] #List to store TPM scaled to the million
  
        for index, row in count_df.iterrows():
            rg = row.iloc[6] #Gene counts
            flg = row['Length']/10**3 #Gene length in kilobases
            TPM_calc = rg/flg #Gene counts normalized by gene length
            TPM_up.append(TPM_calc*1000000) #Normalized gene counts per million
            T += TPM_calc #Add normalized gene counts to the sum
            if check: #If the list of locus tags is empty
                locus_tags.append(row['Geneid']) #Add the locus tag of the gene

        basename = os.path.basename(file) #Get the file name
        TPM_list = [TPM_val/T for TPM_val in TPM_up] #Create a list with the final TPM values by dividing by the sum
        sample = f'I{basename[17:19]}_{basename[20:21]}_{basename[24:27].replace("_", "")}' #Get the name of the sample
        print(f'Calculated TPM for sample {sample}')
        TPM[sample] = TPM_list #Assign the TPMs to the sample

# =============================================================================
# 4. Save TPM calculations to files
# =============================================================================

df = pd.DataFrame.from_dict(TPM, orient='index').transpose() #Convert the TPM dictionary to a dataframe
df['locus_tag'] = locus_tags #Create a new column with locus tags
df = df.set_index('locus_tag') #Convert the new column into an index
df.to_csv(outfile, sep = '\t') #Save the dataframe to the output file as TSV
print('TPM saved to file') 

# =============================================================================
# 5. Calculate mean TPMs and save to file
# =============================================================================

isol_cond = [f'{i}{c}' for i in ['01', '02', '09', '10'] for c in ['F', 'S']] #Loop through isolates and media conditions
mean_df = pd.DataFrame() #Initialize an emoty dataframe
 
for ic in isol_cond: #Loop through isolates+condition
    underscore_ic = ic[:-1] + '_' +  ic[-1:] #Add an underscore between isolate and condition to match dataframe headers
    mean_df[ic] = df.filter(regex=underscore_ic).mean(axis=1) #Apply the mean function on columns with the same isolate+condition in the header
    mean_df.to_csv(mean_outfile, sep = '\t') #Write dataframe to file
print('Mean values saved to files')

# =============================================================================
# 6. Define variable with genes to annotate manually
# =============================================================================

replace_dict = {'AKUH3B104J_00510': 'adhesin_510', 'AKUH3B104J_00520': 'GT2_520',
                'AKUH3B104J_01020': 'adhesin_1020', 'AKUH3B104J_12990': 'GS',
                'AKUH3B104J_13000': 'BrS_13000', 'AKUH3B104J_13010': 'BrS_13010',
                'AKUH3B104J_13020': 'GS-BrS', 'AKUH3B104J_13060': 'wzx',
                'AKUH3B104J_13080': 'wzy', 'AKUH3B104J_13090': 'GT14_13090',
                'AKUH3B104J_13100': 'GT14_13100', 'AKUH3B104J_13110': 'GT2_13110', 
                'AKUH3B104J_13120': 'GT1_13120', 'AKUH3B104J_13130': 'epsE', 
                'AKUH3B104J_13140': 'epsD', 'AKUH3B104J_13150': 'epsC', 
                'AKUH3B104J_13160': 'epsB', 'AKUH3B104J_13170': 'epsA',
                'AKUH3B104J_14310': 'adhesin_14310', 
                'AKUH3B104J_PKUN00010': 'parA', 'AKUH3B104J_PKUN00020': 'repC',
                'AKUH3B104J_PKUN00030': 'kukP', 'AKUH3B104J_PKUN00040': 'kukA',
                'AKUH3B104J_PKUN00050': 'kukC', 'AKUH3B104J_PKUN00060': 'kukT',
                'AKUH3B104J_PKUN00070': 'kukF', 'AKUH3B104J_PKUN00080': 'kukE',
                'AKUH3B104J_PKUN00090': 'kukG_90', 'AKUH3B104J_PKUN00100': 'kukG_100',
                'AKUH3B104J_PKUN00110': 'kukB', 'AKUH3B104J_PKUN00120': 'pKUN_HK',
                'AKUH3B104J_PKUN00130': 'pKUN_HTH', 'AKUH3B104J_PKUN00140': 'pKUN_tnpR',
                'AKUH3B104J_PKUN00190': 'repA'}

# =============================================================================
# 7. Update the annotations
# =============================================================================

print(f'Step 0: Parsing GenBank file: {ref_gbk}...')
ref_dict = parse_gbk(base_gbk, ref_gbk)
print('Step 0: Done!') 

print(f'Parsing comparison {mean_outfile}')
print(f'Step 1/4: Open {mean_outfile}')
gene_names = []
descriptions = []
with open(mean_outfile) as TPM_csv:
    print(f'Step 2/4: Read {mean_outfile} as a dataframe')
    TPM_df = pd.read_csv(TPM_csv, sep = '\t', index_col = 0)
    print('Step 3/4: Loop through dataframe and update the annotations')
    for index, row in TPM_df.iterrows(): #Loop through the dataframe
        if index in ref_dict.keys():
            descriptions.append(ref_dict[index][1])
            if index not in replace_dict.keys():
                gene_names.append(ref_dict[index][0])
            else:
                gene_names.append(replace_dict[index])
        else: descriptions.append('-'), gene_names.append('-')
    
    TPM_df.insert(0, 'Description', descriptions) #Insert annotation column at the front
    TPM_df.insert(1, 'Gene_name', gene_names) #Insert annotation column at the front
            
with open(annot_outfile, 'w') as out_tsv:
    print(f'Step 4/4: Write to output file: {annot_outfile}')
    TPM_df.to_csv(out_tsv, sep = '\t')
print('4/4 steps done!')