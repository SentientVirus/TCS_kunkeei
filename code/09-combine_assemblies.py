#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 15:24:14 2025

This script takes the HGAP chromosome assemblies from NGI and the plasmid 
sequences from Flye assemblies and merges them into single files for each
isolate.

@author: Marina Mota-Merlo
"""

import os
import logging, traceback
from Bio import SeqIO

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
# 1. Define inputs from Snakemake
# =============================================================================

og_strain = snakemake.input.og_strain #Reference genome assembly
NGI_inputs = snakemake.input.NGI #NGI assemblies
Flye_inputs = snakemake.input.Flye #Flye assemblies
outfiles = snakemake.output.assemblies #Combined assemblies
outdir = os.path.dirname(outfiles[0]) #Path to the combined assemblies
mauve_dir = os.path.dirname(snakemake.output.pgvmauve) #Path to the progressiveMauve output

#Create output directory for combined assemblies if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

#Create .fasta files with combined assemblies
print('Loop through files...')
for i in range(0, len(NGI_inputs)): #Loop through isolates using an index
    print(f'Retrieve chromosome from file {NGI_inputs[i]}...')
    chromosome = list(SeqIO.parse(NGI_inputs[i], 'fasta'))[0] #Retrieve the first sequence (chromosome)
    chromosome.id = 'chromosome' #Rename contig representing the chromosome
    chromosome.description = 'chromosome' #Update description
    
    print(f'Retrieve plasmid from file {Flye_inputs[i]}...')
    plasmid = list(SeqIO.parse(Flye_inputs[i], 'fasta'))[1] #Retrieve the second sequence (plasmid)
    
    print(f'Write sequences to {outfiles[i]}...')
    with open(outfiles[i], 'w') as out_fasta: #Open the output file in write mode
        SeqIO.write([chromosome, plasmid], out_fasta, 'fasta') #Write the chromosome and plasmid information
        
#Create a list with all mauve input files separated by spaces
mauve_inputs = f'{og_strain} ' 
for file in outfiles:
    mauve_inputs += f'{file} '
    
#Run progressiveMauve on the combined assemblies
os.system(f'pgv-pmauve {mauve_inputs[:-1]} -o {mauve_dir} >> {snakemake.log[0]} 2>> {snakemake.log[0]}')
    
