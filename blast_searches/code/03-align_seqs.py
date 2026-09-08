#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 14:27:43 2026

Script to run MAFFT L-INS-i to align the retrieved amino acid sequences to
each other.
Environment: pixi_phylo/default.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import subprocess
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
# 1. Set paths to inputs and outputs
# =============================================================================

threads = snakemake.threads #No. of threads to be used by MAFFT

for i in range(len(snakemake.output)):

    in_faa = snakemake.input[i] #Path to the formatted fasta file
    out_faa = snakemake.output[i] #Path to the alignment file
    
# =============================================================================
# 2. Run MAFFT
# =============================================================================
    
    print(f'Running MAFFT on {in_faa} and saving results to {out_faa}... (1/1)')
    #Define the command to run
    command = f'mafft-linsi --thread {threads} {in_faa} > {out_faa} 2> {log}'
    subprocess.run(command, shell = True) #Run MAFFT
    print('Done! (1/1)')