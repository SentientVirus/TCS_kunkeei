#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 14:07:57 2025

Script that modifies the NCBI IDs of the CDS in the GFF file so that it can be
used with emapper2gbk.

@author: Marina Mota-Merlo
"""

import os
import logging, traceback

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
# 1. Set paths to inputs and outputs
# =============================================================================

input_gff = snakemake.input[0]
output_gff = snakemake.output[0]

# =============================================================================
# 2. Read the input and create the output
# =============================================================================

print(f'Open input: {input_gff} and output: {output_gff}')
with open(input_gff) as in_gff, open(output_gff, 'w') as out_gff: #Open the input and output files
    print('Looping through input...')
    i = 1
    for line in in_gff.readlines(): #Loop through the lines in the input file
        print(f'Reading line {i}...')
        line = line.replace('cds-', '') #In each line, remove the substring 'cds-'
        out_gff.write(line) #Write the modified line to the output file
        i += 1
    print('Done!')
                

    