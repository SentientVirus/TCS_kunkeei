#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 14:39:03 2026

Script to run IQtree on the MAFFT L-INS-i alignment to generate a phylogeny of
the RR-TF.
Environment: pixi_phylo/default.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import subprocess
import logging, traceback
import sys
import os

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0]

with open(log, 'w') as handle: #Open the log file in write mode
    handle.write('') #Overwrite the file

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
# 1. Set paths to inputs and outputs
# =============================================================================

outpath = os.path.dirname(snakemake.output[0]) #Directory to store the outputs
threads = snakemake.threads #No. of threads to be used by IQtree

if not os.path.exists(outpath): #If the output directory doesn't exist
    os.makedirs(outpath) #Create it


for i in range(len(snakemake.input)):

    in_aln = snakemake.input[i] #Path to the formatted fasta file

# =============================================================================
# 2. Run IQtree
# =============================================================================
    
    print(f'Running IQtree on {in_aln} and saving results to {outpath}... (1/1)')
    #Define the command to run
    command = f'iqtree -nt AUTO -ntmax {threads} -s {in_aln} -st AA -msub nuclear -bb 1000 -bnni >> {log}'
    subprocess.run(command, shell = True) #Run IQtree
    subprocess.run(f'mv {in_aln}.* {outpath}', shell = True) #Move IQtree results to the desired folder
    print('Done! (1/1)')