#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:15:16 2026

Script to divide the file with all the proteins from strain H3B1-04J into
several files, each with a maximum of 100 proteins, to run the InterProScan
online server.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 1. Import required modules
# =============================================================================

import os
from Bio import SeqIO

# =============================================================================
# 2. Define paths to inputs and outputs
# =============================================================================

gbk = os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff' #GenBank input (to retrieve locus tags)
faa = gbk.replace('genomic', 'protein').replace('gbff', 'faa') #Path to the file with the protein sequences in FASTA format
workdir = os.path.expanduser('~') + '/mucoid_project/snpseq00064' #Working directory
outdir = f'{workdir}/data/interpro' #Directory to store outputs

#Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

# =============================================================================
# 2. Create a dictionary with protein ID: locus tag pairs
# =============================================================================

loctag_dict = {} #Empty dictionary
with open(gbk) as handle: #Open GenBank file
    for record in SeqIO.parse(handle, 'genbank'): #Loop through the records (genomic elements) in the GenBank
        record_dict = {feature.qualifiers['protein_id'][0]: feature.qualifiers['locus_tag'][0] for feature in record.features if ('locus_tag' in feature.qualifiers.keys()) and ('protein_id' in feature.qualifiers.keys())} #Loop through the CDS and retrieve the protein ID: locus tag pairs
        loctag_dict = record_dict | loctag_dict #Merge the dictionary with the previous dictionary, so that the proteins in both genomic elements are included

# =============================================================================
# 3. Divide the original FASTA file into smaller files with <=100 sequences        
# =============================================================================

i = 1 #Index of the protein sequence
n = 1 #Index indicating the file where the sequence will be saved
to_write = [] #Empty list

with open(faa) as handle: #Open the fasta file with protein sequences
    for record in SeqIO.parse(handle, 'fasta'): #Loop through the protein records in the file
        record.description = record.description.replace(f'{record.id} ', '') #Remove the protein ID from the description
        record.description = record.description.replace(f'{loctag_dict[record.id]} ', '') #Remove the locus tag from the description
        record.id += f'/{loctag_dict[record.id]}' #Add the locus tag to the record ID (protein ID)
        
        if i > 100: #If there are more than 100 sequences in a file
            n += 1 #Increase the file index
            i = 1 #Reset i
            
        outfile = f'{outdir}/H3B1-04J_proteins{n}.faa' #Set the output file
        
        if i == 1: #If it is the same protein in the file
            with open(outfile, 'w') as faa_out: #Overwrite the file
                SeqIO.write(record, faa_out, 'fasta') #Write the protein record
        else: #Otherwise
            with open(outfile, 'a') as faa_out: #Append to the file
                SeqIO.write(record, faa_out, 'fasta') #Write the protein record
        i += 1 #Increase the index of the next protein sequence