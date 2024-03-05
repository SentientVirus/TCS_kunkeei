# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 18:24:07 2021

Script to Blast all the genes, annotated by Prokka, from the newly-sequenced
genomes to the reference genome.

@author: Marina Mota-Merlo
"""
import os
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
# Defining inputs
# =============================================================================

files = snakemake.input.fnas
subject = snakemake.input.og_strain
outpath = snakemake.params.outpath
threads = snakemake.threads

# =============================================================================
# Running Blast       
# =============================================================================

for i in range(len(files)):
    S = f'blastn -query {files[i]} -out {outpath}/{files[i][:-4].split("/")[2]}.tab \
        -subject {subject} -outfmt 7'
    os.system(S)
    
