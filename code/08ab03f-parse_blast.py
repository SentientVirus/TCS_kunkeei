#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:43:35 2024

This script will go through the all-vs-all BLAST results and return a list of
the genes that don't have any identical match in the reference strain.

@author: Marina Mota-Merlo
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

    logger.error(''.join(['Uncaught exception: ',
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# 1. Defining inputs 
# =============================================================================

infiles = snakemake.input #Path to input files
outfile = snakemake.output[0] #Path to output files

# =============================================================================
# 2. Function to retrieve unique locus tag in each sample
# =============================================================================

def get_unique_seqs(file, outdict = {}, id_dict = {}):
    '''Function to retrieve genes that are different between the reference
    assembly and the re-sequenced strains
    
    Inputs:
        - file: Output file from BLAST (outfmt 6).
        - outdict: Dictionary to store the locus tags of genes that differ to
                   the reference genome.
        - id_dict: Dictionary that assigns the reference protein ID to each 
                   locus tag.'''
        
    bool_dict = {} #Initialize empty dictionary
    print(f'Processing {file}...')
    filename = os.path.basename(file).split('.')[0] #Retrieve the file name
    
    with open(file) as infile: #Open input BLAST
        df = pd.read_csv(infile, sep = '\t', header = None) #Read file as dataframe
        for index, row in df.iterrows(): #Loop through rows in the dataframe
            if row[0] not in bool_dict.keys(): #If the locus tag of the gene is not in the boolean dictionary
                bool_dict[row[0]] = True #Add it to the dictionary with a True value
                id_dict[(filename, row[0])] = row[1].split('_')[2] #Assign the protein ID from the reference genome to the locus tag
            if row[2] == 100 and row[3] == row[9]: #If the match has 100% identity and the length of the genes is identical
                bool_dict[row[0]] = False #Set the boolean value to False
            
    true_list = [k for k in bool_dict.keys() if bool_dict[k] == True] #Retrieve locus tags of genes that differ to the reference
    outdict[filename] = true_list #Assign the list of locus tags to the file name

# =============================================================================
# 3. Function implementation
# =============================================================================
true_dict = {} #Dictionary to store locus tags of genes that differ
loc_dict = {} #Dictionary to map the reference protein IDs to the locus tags
[get_unique_seqs(file, true_dict, loc_dict) for file in infiles]; #Implement the function on all input files

# =============================================================================
# 4. Saving results to file
# =============================================================================
with open(outfile, 'w') as txt: #Open output file in write mode
    txt.write('sample_no\tlocus_tag\treference_locus\n') #Write headers
    [txt.write(f'{sample}\t{locus}\t{loc_dict[(sample, locus)]}\n') for (sample, loci) in true_dict.items() for locus in loci] #Write information to file
    
# =============================================================================
# 5. Retrieving the protein ID of candidate loci
# =============================================================================

with open(outfile) as handle: #Open previous output file
    dif = pd.read_csv(handle, sep = '\t') #Read it as dataframe
    
samples = list(set(dif['sample_no'])) #List of samples (isolates)

#Create a dictionary to assign protein IDs to each sample
prot_dict = {sample: list(dif[dif['sample_no'] == sample]['reference_locus']) for sample in samples}

#Get protein IDs that are found only in the two first samples (mucoid)
possible_dif = [gene for gene in prot_dict[samples[0]] if gene in prot_dict[samples[1]] and gene not in prot_dict[samples[2]] + prot_dict[samples[3]]]

print('List of mutated genes that could explain the deletion:')
[print(gene) for gene in possible_dif]
