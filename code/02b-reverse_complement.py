#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:26:50 2023

Script to reverse the strand of the genome of isolate 10, which has the
reverse strand as forward strand in the raw assembly.

@author: Marina Mota-Merlo
"""

from Bio import SeqIO
import logging, traceback

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

filename = str(snakemake.input) #Input assembly
outfile = str(snakemake.output) #Output assembly with reversed contig
record_list = [] #List of records

# =============================================================================
# Reverse strand and save to file
# =============================================================================

with open(filename) as reverse_stranded: #Open input file
    records = SeqIO.parse(reverse_stranded, 'fasta') #Read records in file as fasta
    for record in records: #Loop through records in file
        record.seq = record.seq.reverse_complement() #Convert the sequences to their reverse complement
        record_list.append(record) #Append modified records to list

with open(outfile, 'w') as output: #Open output file in write mode
    SeqIO.write(record_list, output, 'fasta') #Write modified records to file
