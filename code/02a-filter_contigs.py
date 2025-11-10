#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:26:50 2023

Script to filter out contigs with low coverage, which are unlikely to be true
genomic elements.

@author: Marina Mota-Merlo
"""

import os
from Bio import SeqIO
import logging, traceback
import pandas as pd

# =============================================================================
# Logging
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
# Define inputs and outputs
# =============================================================================

infiles = snakemake.input.stats #Files with assembly statistics
assemblies =  snakemake.input.assembly #Flye assemblies
outfiles = snakemake.output #Output assembly files without the low-coverage contigs

# =============================================================================
# Filter out contigs and save to file
# =============================================================================

for i in range(len(infiles)): #Loop through input files
    to_exclude = [] #Create an empty list to add contigs to filter out
    with open(infiles[i]) as handle: #Open input file
        print('Filename: ', infiles[i]) #Print filename
        df = pd.read_csv(handle, sep = '\t') #Read file as dataframe
        for index, contig in df.iterrows(): #Loop through rows in the dataframe
            contig_name = contig['#seq_name'] #Retrieve contig name
            print('Contig name: ', contig_name) #Print contig name
            if contig['cov.'] < 30 or int(contig['length']) < 18000: #If the coverage of the contig is below 30
                to_exclude.append(contig_name) #Append it to the list of contigs to exclude
        print('Contigs to exclude: ', to_exclude) #Print the list
        
    with open(assemblies[i]) as assembly: #Open the corresponding assembly
        records = SeqIO.parse(assembly, 'fasta') #Read the fasta records in the file
        to_write = [record for record in records if record.id not in to_exclude] #Create a list of records to keep
        for record in to_write:
            if len(record.seq) > 20000:
                record.id = 'chromosome'
                record.description = 'chromosome'
            else: 
                record.id = 'pKUN'
                record.description = 'pKUN'
                
    with open(outfiles[i], 'w') as outfile: #Open output file
        SeqIO.write(to_write, outfile, 'fasta') #Write records to file
