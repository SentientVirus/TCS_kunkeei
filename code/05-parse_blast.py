#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:43:35 2024

This script will go through the all-vs-all BLAST results and return a list of
the genes that don't have any identical match in the reference strain.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Importing packages
# =============================================================================
import os
import logging, traceback
import pandas as pd

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

    logger.error(''.join(['Uncaught exception: ',
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# 1. Defining inputs 
# =============================================================================
infiles = snakemake.input
outfile = snakemake.output

# =============================================================================
# 2. Function to retrieve unique locus tag in each sample
# =============================================================================

def get_unique_seqs(file, outdict = {}):
    bool_dict = {}
    with open(file) as infile:
        df = pd.read_csv(infile, sep = '\t', header = None)
        for index, row in df.iterrows():
            if row[0] not in bool_dict.keys():
                bool_dict[row[0]] = True
            if row[3] == row[9]:
                bool_dict[row[0]] = False
            
    true_list = [k for k in bool_dict.keys() if bool_dict[k] == True]
    filename = os.path.basename(file).split('.')[0]
    outdict[filename] = true_list

# =============================================================================
# 3. Function implementation
# =============================================================================
true_dict = {}
[get_unique_seqs(file, true_dict) for file in infiles];

# =============================================================================
# 4. Saving results to file
# =============================================================================
with open(outfile, 'w') as txt:
    txt.write('sample_no\tlocus_tag\n')
    [txt.write(f'{sample}\t{locus}\n') for (sample, loci) in true_dict.items() for locus in loci]
