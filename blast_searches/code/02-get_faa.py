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
import logging
import sys

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0]

#Overwrite the log file every time
with open(log, 'w') as logfile:
    logfile.write('')
            
#Redirect stdout and stderr to log file
sys.stdout = open(log, 'a')
sys.stderr = open(log, 'a')

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

# =============================================================================
# 1. Set path to output
# =============================================================================

print('test')
for i in range(len(snakemake.output.faa)):

    outfile = snakemake.output.faa[i] #Output fasta file with the protein sequences
    finalfile = snakemake.output.formatted[i] #Same, but the annotation of the protein sequences is shortened to later include it in the tree
    descriptions = snakemake.input[i] #File with the BLASTp hit descriptions

    print(f'Processing {descriptions} to generate {outfile} and {finalfile}...')
    
    # =============================================================================
    # 2. Define the sequences to be downloaded
    # =============================================================================
    
    if 'HK' in outfile:
        #List of hits to be retrieved after reviewing the BLASTp results
        acc_list = ['WP_353317292.1', 'WP_120784741.1', 'WP_203619427.1', 
                    'WP_056966301.1', 'WP_353485676.1', 'WP_220728880.1', 
                    'WP_220729921.1', 'CAI2564375.1', #Apilactobacillus
                    
                    'WP_057974193.1', 'WP_039145019.1', 'WP_128517339.1', 
                    'WP_211761223.1', 'WP_057771063.1', 'WP_413526435.1', 
                    'WP_003552833.1', 'WP_045807370.1', 'WP_270333345.1', 
                    'WP_229267429.1', 'WP_224288406.1', 'WP_146990902.1',
                    'WP_150204492.1', 'WP_457313907.1', 'WP_057818988.1', 
                    'WP_155431510.1', #Other Lactobacillaceae
                    
                    'WP_101703914.1', 'WP_069698772.1', 'WP_193991767.1',
                    #Other Lactobacillales (first is Carnobacteriaceae, rest are Enterococcaceae)
                    'WP_014525065.1', 'WP_001253704.1', 'WP_003917137.1', 
                    #Model CroS, EnvZ and MtrB
            ]
    elif 'RR-TF' in outfile:
        acc_list = ['WP_125712160.1', 'WP_353317293.1', 'WP_120784742.1',
                    'WP_203619426.1', 'WP_054658628.1', 'WP_220728881.1', 
                    'CAI2564285.1', #Apilactobacillus
                    
                    'WP_125761468.1', 'WP_089939307.1', #'WP_446183611.1',
                    'WP_069698773.1', 'WP_349641026.1', 'WP_039145020.1',
                    'WP_010690715.1', 'WP_407884581.1', 'WP_056974007.1',
                    'WP_128517338.1', 'WP_211761224.1', 'WP_057771000.1',
                    'WP_413526436.1', 'WP_003552832.1', 'WP_048592711.1',
                    'WP_045807369.1', 'WP_270333343.1', 'WP_229267428.1',
                    'WP_224288405.1', 'WP_146990904.1', 'WP_339959950.1', 
                    'WP_459055447.1', 'WP_150204493.1', 'WP_457313909.1',
                    'WP_168721664.1', 'WP_057818986.1', 'WP_155431509.1',
                    'WP_063084080.1', #Other Lactobacillaceae
                    
                    'WP_101703913.1', 'WP_013774648.1', 'WP_028790024.1',
                    'WP_077275945.1', #Other Lactobacillales (first is Carnobacteriaceae, rest are Enterococcaceae)
                    'WP_002355963.1', 'WP_001157757.1', 'WP_003899985.1'  
                    #Model CroR, OmpR and MtrA
            ]
    
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
            elif 'CroS' in record.description:
                gtype = 'CroS'
            else: gtype = ''
    
            if record.id == 'WP_014525065.1':
                genus = 'Enterococcus'
                epithet = 'model '
            elif record.id == 'WP_002355963.1':
                genus = 'Bacterial'
                epithet = 'model '
            elif record.id == 'WP_001253704.1':
                genus = 'Escherichia'
                epithet = 'coli '
                gtype = 'EnvZ'
            elif record.id == 'WP_003917137.1':
                genus = 'Mycobacterium'
                epithet = 'model '
                gtype = 'MtrB'
            elif record.id == 'WP_001157757.1':
                genus = 'Escherichia'
                epithet = 'coli '
                gtype = 'OmpR'
            elif record.id == 'WP_003899985.1':
                genus = 'Mycobacterium'
                epithet = 'model '
                gtype = 'MtrA'
            
            elif len(species.split(' ')) > 1: #If the species name contains more than one word
                epithet = species.split(' ')[1] #Retrieve the species epithet
            else: genus = genus[:-1] #If the genus is the last word, remove the empty space at the end
            for i in range(2, len(species.split(' '))): #Loop through the words after the species
                epithet += f'_{species.split(" ")[i]}' #Add them after the species
            epithet = epithet[:-1] #Remove the character at the end
            
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
