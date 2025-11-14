#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 14:02:19 2025

Script that takes the reference GenBank file and the re-annotated one and adds
the locus tags to the re-annotated one.

@author: Marina Mota-Merlo
"""

import os
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
# 1. Define paths to inputs and outputs
# =============================================================================

infile = snakemake.input.emapper #GenBank file from eggnogmapper
inref = snakemake.input.NCBI #Reference GenBank from NCBI
outfile = snakemake.output[0] #Output file

# =============================================================================
# 2. Function to create a dictionary mapping locus tags to protein ids
# =============================================================================

def retrieve_loctags(ref: str):
    prot2loctag = {} #Create empty dictionary
    with open(ref) as gbk: #Open reference GenBank
        #Create output dictionary using dictionary comprehension
        prot2loctag = {feature.qualifiers['protein_id'][0]: feature.qualifiers['locus_tag'][0] for contig in SeqIO.parse(gbk, 'genbank') for feature in contig.features if 'locus_tag' in feature.qualifiers.keys() and 'protein_id' in feature.qualifiers.keys()}
    return prot2loctag #Return output dictionary

# =============================================================================
# 3. Run the code
# =============================================================================
print(f'Retrieving locus tags from {inref}...')
id2loctag = retrieve_loctags(inref) #Apply the function

print('Done!')

print(f'Reading input file {infile}...')
with open(infile) as gbk: #Open input GenBank from eggNOG-mapper
    lines = '' #Create empty string
    for line in gbk.readlines(): #Loop through lines in the GenBank
        if 'locus_tag' in line: #If the lines contain "locus_tag"
            prot_id = line.replace('/locus_tag=', '').strip() #Retrieve the protein ID from the line (CAI...)
            pid = prot_id.replace('"', '') #Remove extra quotation marks from the string
            line = line.replace(pid, id2loctag[pid]) #Replace the protein ID with the locus tag
        elif 'product' in line: #If the description of the gene is in the line
            line += ' '*21 + f'/protein_id={prot_id}\n' #Add another line under it with the protein ID
        lines += line #Add it to the string created at the start

print('Done!')        

print(f'Creating output file {outfile}...')
with open(outfile, 'w') as new_gbk: #Open output file
    new_gbk.write(lines) #Write the string containing the file contents to the file

print('Done!')