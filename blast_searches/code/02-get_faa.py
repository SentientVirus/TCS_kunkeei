#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 17:51:55 2025

Code to retrieve protein sequences from other A. kunkeei species from NCBI,
add the species name to the accession, and save the sequences to a single
fasta file.
Environment: pixi_phylo/default.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import subprocess
import os
from Bio import SeqIO
import pandas as pd
import logging, traceback
import sys

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0]

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    '''Function to handle exceptions'''
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger = logging.getLogger() #Create a logger

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ])) #Format to show exceptions

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Set path to output
# =============================================================================

outfile = snakemake.output.faa #Output fasta file with the protein sequences
finalfile = snakemake.output.formatted #Same, but the annotation of the protein sequences is shortened to later include it in the tree
descriptions = snakemake.input[0] #File with the BLASTp hit descriptions

# =============================================================================
# 2. Define the sequences to be downloaded
# =============================================================================

#List of hits to be retrieved after reviewing the BLASTp results
acc_list = ['WP_353317293.1', 'WP_120784742.1', 'WP_054658628.1', #Apilactobacillus
            'WP_220728881.1', 'CAI2564285.1',
            'WP_035452186.1', 'WP_125761468.1','WP_089939307.1', #Other Lactobacillaceae
            'WP_349641026.1', 'WP_395391563.1', 'WP_010690715.1', 
            'WP_407884581.1', 'WP_056974007.1', 'WP_123156317.1', 
            'WP_280135435.1', 'WP_003627353.1', 'WP_125649329.1', 
            'WP_413526436.1', 'WP_218711211.1', 'WP_048592711.1', 
            'WP_462926425.1', 'WP_270333343.1', 'WP_465068414.1', 
            'WP_224288405.1', 'WP_219480889.1', 'WP_339959950.1', 
            'WP_459055447.1', 'WP_150204493.1', 'WP_199877929.1', 
            'WP_168721664.1', 'WP_057818986.1', 'WP_155431509.1', 
            'WP_028790024.1', 'WP_063084080.1',
            'WP_034561863.1', 'MGX7713849.1', 'WP_013774648.1', #Other non-Lactobacillaceae
            'WP_077275945.1', 'WP_048940986.1']

if not os.path.exists(os.path.dirname(outfile)): #If the path to the output file doesn't exist
    os.makedirs(os.path.dirname(outfile)) #Create it
    
with open(outfile, 'w') as ncbi_out: #Open the output file
    ncbi_out.write('') #Overwrite file contents
    
desc_df = pd.read_csv(descriptions, sep = ',') #Read hit descriptions as a dataframe

# =============================================================================
# 3. Create a string with the command to be run
# =============================================================================

list_command = '' #Create an empty string
for acc in acc_list: #Loop through accessions to retrieve
    list_command += f'{acc} ' #Add them to the string
list_command = list_command.strip() #Remove the space at the end of the string

print(f'Records to be downloaded: {list_command}')

#Prepare the complete command to download NCBI records
full_command = f'ncbi-acc-download --molecule protein --out {outfile} {list_command}'

# =============================================================================
# 4. Run the command
# =============================================================================

print('Run ncbi-acc-download... (1/3)')
subprocess.run(full_command, shell = True)
print('Done! (1/3)')

# =============================================================================
# 5. Read the output file and modify the fasta headers
# =============================================================================
    
print('Loop through the downloaded records to format them... (2/3)')
records = [] #Create an empty list
with open(outfile) as handle: #Open the output file
    for record in SeqIO.parse(handle, 'fasta'): #Loop through the records
        print(f'Reading {record.id}...')
        record.description = record.description.replace('[[', '[').replace('] ', ' ') #Retrieve the strings with the species name
        species = record.description.split('[')[1] #Get the genus & species name
        genus = species.split(' ')[0] #From there, retrieve the genus
        if 'CroR' in record.description: #If the record description contains CroR
            gtype = 'CroR' #Assign it as the gene type
        else: gtype = ''
        if len(species.split(' ')) > 1: #If the species name contains more than one word
            epithet = species.split(' ')[1][:-1] #Retrieve the species epithet
        if len(species.split(' ')) > 2 and epithet == 'sp': #If the species epithet is sp
            epithet += f'_{species.split(" ")[2]}_{species.split(" ")[3][:-1]}' #Add an extra word to the species name
        elif len(species.split(' ')) > 2:
            epithet += f'i_{species.split(" ")[2][:-1]}'
        elif len(species.split(' ')) == 1: #If the species is only one word
            epithet = '' #Leave the species epithet empty
            genus = genus[:-1] #Retrieve the genus
        new_description = f'{record.id}_{genus}_{epithet}_{gtype}' #Define new record description
        if new_description.endswith('_'): #If the description ends with an underscore
            new_description = new_description[:-1] #Remove it
        record.description = new_description #Assign the new description to the record
        record.id = new_description #Use it as the record id too
        print(f'New record description: {new_description}') #Print description
        records.append(record) #Append the record to the list of records
print('Done! (2/3)')
        
print(f'Writing formatted records to {finalfile}... (3/3)')
with open(finalfile, 'w') as handle: #Open the output file
    SeqIO.write(records, handle, 'fasta') #Write the records to the file
print('All done!')
