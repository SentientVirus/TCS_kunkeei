#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 11:05:32 2023

Script to run progressive Mauve to compare the genomes of each isolate to the
previous one, including the reference isolate.

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
# Defining input files
# =============================================================================

infiles = snakemake.input['og_strain'] + ' '
for infile in snakemake.input['new_seqs']:
    infiles += f'{infile} '

outdir = os.path.dirname(snakemake.output[0])
    
# =============================================================================
# Running progressive Mauve
# =============================================================================
        
os.system(f'pgv-pmauve --seq_files {infiles}\
            -o {outdir} --tick_style bar >> {snakemake.log[0]}')


