#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 17:45:37 2025

This script reads the output files from DESeq2 and adds annotations from the
GenBank file with refined annotations from the genome re-sequencing project.

@author: Marina Mota-Merlo
"""

import os
import pandas as pd
from Bio import SeqIO
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
# 1. Define functions to be used in the script
# =============================================================================

def parse_gbk(gbk: str):
    """
    Function that takes the path to a GenBank file and parses its contents to
    create a dictionary that maps gene names to locus tags.

    Parameters
    ----------
    gbk : str
        Path to the input GenBank file.

    Returns
    -------
    ref_dict : dict
        A dictionary where locus tags are the keys and gene name annotations
        are the values.

    """
    
    ref_dict = {} #Create output dictionary
    with open(ref_gbk) as gbk: #Open input GenBank
        for contig in SeqIO.parse(gbk, 'genbank'): #Loop through records (contigs) in the GenBank
            for CDS in contig.features: #Loop through features
                if 'locus_tag' in CDS.qualifiers.keys(): #If the feature has a locus tag
                    loctag = CDS.qualifiers['locus_tag'][0] #Retrieve it
                    if 'gene' in CDS.qualifiers.keys(): #If the feature has the "gene" property (gene name)
                        gname = CDS.qualifiers['gene'][0] #Retrieve it
                    else: gname = '-' #Else, set "-" as gene name
                    
                    ref_dict[loctag] = gname #Assign the gene name to the locus tag in the output dictionary
                    
    return ref_dict

# =============================================================================
# 2. Define inputs and outputs from Snakemake
# =============================================================================

ref_gbk = snakemake.input.ref_gbk #os.path.expanduser('~') + '/mucoid_project/ugc00027/results/annotations/emapper2gbk/reference_loctag.gbk'
#indir = os.path.expanduser('~') + '/mucoid_project/snpseq00064/results/DE'
infiles = snakemake.input.DE_annot #[f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('_annotated.tsv')]
outfiles = snakemake.output #[file.replace('_annotated.tsv', '_improved_annot.tsv') for file in infiles]

print(f'Reference GenBank: {ref_gbk}')

# =============================================================================
# 3. Define variable with genes to annotate manually
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
                'AKUH3B104J_PKUN00040': 'kukA', 'AKUH3B104J_PKUN00110': 'nisB'}

# =============================================================================
# 4. Update the annotations
# =============================================================================

print(f'Step 0: Parsing GenBank file: {ref_gbk}...')
ref_dict = parse_gbk(ref_gbk)
print('Step 0: Done!') 


for i in range(0, len(infiles)): #Loop through input files
    comparison = os.path.basename(infiles[i]) #Retrieve filename of the input files
    comparison = comparison.replace('_lfc1', '').replace('_annotated.tsv', '').replace('_', ' ') #Retrieve comparison name
    print(f'Parsing comparison {comparison}')
    print(f'Step 1/4: Open {infiles[i]}')
    with open(infiles[i]) as DE_csv:
        print(f'Step 2/4: Read {infiles[i]} as a dataframe')
        DE_df = pd.read_csv(DE_csv, sep = '\t', index_col = 0)
        print('Step 3/4: Loop through dataframe and update the annotations')
        for index, row in DE_df.iterrows(): #Loop through the dataframe
            gname = ref_dict[index]
            DE_df.at[index, 'gene_names'] = gname
            if index in replace_dict:
                DE_df.at[index, 'gene_names'] = replace_dict[index]
    with open(outfiles[i], 'w') as out_tsv:
        print(f'Step 4/4: Write to output file: {outfiles[i]}')
        DE_df.to_csv(out_tsv, sep = '\t')
    print('4/4 steps done!')