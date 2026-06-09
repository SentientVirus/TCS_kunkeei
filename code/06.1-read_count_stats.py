#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 19:05:58 2026

Script to parse the .summary files outputted by featureCounts, retrieve general
statistics for each sample and save to file.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, sys
import logging
import pandas as pd

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0] #Path to the log file
logdir = os.path.dirname(log) #Path to the log directory

#Create the log directory if it doesn't exist
if not os.path.exists(logdir):
    os.makedirs(logdir)
    
#Overwrite log file
with open(log, 'w') as logfile:
    logfile.write('')
            
#Redirect stdout and stderr to log file
sys.stdout = open(log, 'a')
sys.stderr = open(log, 'a')

#Logging configuration
logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================

#List of input files
infiles = sorted(snakemake.input)

#Path to output file
outfile = snakemake.output[0]

logging.info(f'Create output directory: {os.path.dirname(outfile)}')
#Create output directory if it doesn't exist
if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))
    
logging.info(f'Create output file: {os.path.basename(outfile)}')
#Overwrite the output file
with open(outfile, 'w') as handle: #Open the file in write mode
    #Write headers to file
    handle.write('Sample\tTotal reads (counts)\tReads mapped (counts)\t')
    handle.write('Reads mapped (%)\tAssigned (counts)\t')
    handle.write('Assigned (% of mapped)\tAssigned (% of total)\n')

# =============================================================================
# 2. Loop through the files and generate the statistics
# =============================================================================

logging.info(f'Loop through the {len(infiles)} input files...')
for file in infiles: #Loop through input files
    
    filename = os.path.basename(file) #Retrieve file name
    
    #Retrieve sample name
    sample = filename.split('_')[0].replace('VF-3336-H3B1-04J-', '')
    
    logging.info(f'Processing sample {sample}...')
    
    #Read input count summary file as a dataframe
    with open(file) as handle:
        count_info = pd.read_csv(handle, sep = '\t')
        
    #Rename the dataframe columns
    count_info.columns = ['Status', 'Count']
    
    #Retrieve the desired statistics
    total_reads = count_info['Count'].sum() #Sum of all reads
    mapped_reads = count_info.drop(1).sum().values[1] #Sum of all reads except unmapped
    assigned_reads = count_info.iloc[0]['Count'] #No. of assigned reads
    
    perc_mapped = mapped_reads/total_reads*100 #% of reads mapped
    mapped_assigned = assigned_reads/mapped_reads*100 #% of mapped reads that are assigned to genes
    perc_assigned = assigned_reads/total_reads*100 #% of all reads that are assigned to genes
    
    #Open the output file and write the information to it
    with open(outfile, 'a') as handle:
        handle.write(f'{sample}\t{total_reads}\t{mapped_reads}\t')
        handle.write(f'{perc_mapped:.2f}\t{assigned_reads}\t')
        handle.write(f'{mapped_assigned:.2f}\t{perc_assigned:.2f}\n')
        
    logging.info(f'Information from {sample} written to file!')
    