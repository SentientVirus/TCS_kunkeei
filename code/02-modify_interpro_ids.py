#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 12:35:32 2025.
Modified from a script created by Marina Mota-Merlo on Thu Apr 13 11:45:45 2023

Script to change the column with protein IDs in the Interproscan output to 
a column with locus tags, so that these files can be navigated using the same
scripts as for the A. kunkeei annotations made by Dyrhage et al. (2022).

@author: Marina Mota-Merlo
"""

# =============================================================================
# Script to change protein IDs in Interproscan annotations to locus tags
# =============================================================================

import os
import pandas as pd
from Bio import SeqIO
import multiprocessing
from functools import partial
import time
import logging, sys

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/adhesins' #Path to the working directory

logdir = f'{workdir}/logs' #Path to the log directory

#Create the log directory if it doesn't exist
if not os.path.exists(logdir):
    os.makedirs(logdir)
    
#Path to the log file
log = f'{logdir}/02-modify_interpro_ids.log'
with open(log, 'w') as logfile: #Overwrite log file
    logfile.write('')
            
#Redirect stdout and stderr to log file
sys.stdout = open(log, 'a')
sys.stderr = open(log, 'a')

#Logging configuration
logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

start_time = time.time() #Get starting time

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================

outpath = f'{workdir}/interproscan/locus_tags' #Output path of InterProScan annotations with locus tags
path = os.path.expanduser('~') + '/Akunkeei_files/gbff' #Input path to GenBank files to retrieve the locus tags from
new_strain_path = path.replace('gbff', 'new_genomes/gbff') #Input path to the GenBank files for "new" strains
prot_loctag_out = f'{workdir}/metadata/prot_id_loctag.tsv' #Path to output metadata file
#Retrieve all the input GBFF files
infiles = [f'{path}/{file}' for file in os.listdir(path) if not os.path.isdir(f'{path}/{file}')] + [f'{new_strain_path}/{file}' for file in os.listdir(new_strain_path)]

no_cores = 12 #Set the no. of cores to run the analysis

#Create output directories if they don't exist
if not os.path.exists(outpath):
    os.makedirs(outpath)
    
if not os.path.exists(os.path.dirname(prot_loctag_out)):
    os.makedirs(os.path.dirname(prot_loctag_out))
    
# =============================================================================
# 2. Function definition
# =============================================================================

def process_infiles(file):
    '''Function to process a GenBank file:
        - Input: 
            file (str): Path to the input file.
        - Output:
            prot2loctag (dict): Dictionary which maps a locus tag to each 
            protein ID.
    '''
    prot2loctag = {} #Create a dictionary to store the locus tag of each protein
    with open(file) as genb: #Open GenBank file
        info = SeqIO.parse(genb, 'genbank') #Parse the file
        strain = os.path.basename(file).split('_')[0] #Retrieve the strain from the file name
        for record in info: #Loop through records (contigs) in the file
            for feature in record.features: #Loop through features (mostly genes) in each record
                tag = feature.qualifiers.get('locus_tag') #Retrieve the locus tag
                prot_id = feature.qualifiers.get('protein_id') #Retrieve the protein ID
                if tag and prot_id: #If there is a locus tag and a protein ID
                    tag = tag[0] #Retrieve the locus tag as string
                    prot_id = prot_id[0] #Retrieve the protein ID as string
                    prot2loctag[prot_id] = tag #Assign the locus tag to the protein ID in the dictionary
    with open(f'{workdir}/interproscan/{strain}.tsv') as df_file: #Open the TSV file with InterProScan results
        df = pd.read_csv(df_file, sep = '\t', header = None) #Read the file as a dataframe
        df2 = df.replace({0: prot2loctag}) #Replace the protein IDs with locus tag
        #Write the modified dataframe to a separate file
        df2.to_csv(f'{outpath}/{strain}.tsv', sep = '\t', header = False, index = False)
    return prot2loctag

# =============================================================================
# 3. Run the function with multithreading
# =============================================================================

logging.info('Process input files...')

#Multithreading
if __name__ == '__main__': 
    pool = multiprocessing.Pool(processes = no_cores) #Create a multiprocessing pool
    outputs = pool.map(partial(process_infiles), infiles) #Run the function
    pool.terminate() #Terminate the pool

logging.info('Write locus tags vs protein IDs...')
with open(prot_loctag_out, 'w') as handle: #Open the output tab file
    handle.write('Locus_tag\tProtein_id\n') #Write headers
    #Write all the locus tag/protein ID pairs
    [{handle.write(f'{value}\t{key}\n') for (key, value) in output.items()} for output in outputs]
    
end_time = time.time() - start_time #Calculate total running time
logging.info(f'This script took {end_time:2f} with {no_cores} processes')
