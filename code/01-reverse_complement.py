#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:26:50 2023

@author: marina
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

filename = snakemake.input
outfile = snakemake.output
record_list = []

# =============================================================================
# Reverse strand and save to file
# =============================================================================

with open(filename) as reverse_stranded:
    records = SeqIO.parse(filename, 'fasta')
    for record in records:
        record.seq = record.seq.reverse_complement()
        record_list.append(record)

with open(outfile, 'w') as output:
    SeqIO.write(record_list, output, 'fasta')
