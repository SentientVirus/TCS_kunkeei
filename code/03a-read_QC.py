#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 17:10:10 2026

Script to run Fastplong on the PacBio reads.

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Import required modules
# =============================================================================

import os, sys
import subprocess
import logging, traceback

# =============================================================================
# 1. Logging
# =============================================================================

log = snakemake.log[0]

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    # Create a logger
    logger = logging.getLogger()

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')


# =============================================================================
# 2. Define inputs and outputs
# =============================================================================

infiles = snakemake.input #List with full paths to input read files
print(infiles)

outfiles = snakemake.output.html
json = snakemake.output.json
reads = snakemake.output.reads
outdir = os.path.dirname(outfiles[0]) #Path to store the FastQC outputs

cores = snakemake.threads #Number of cores

# =============================================================================
# 4. Run FastQC and MultiQC
# =============================================================================

for i in range(len(infiles)):
    
    print(f'Input: {infiles[i]}')
    print(f'Trimmed reads: {reads[i]}')
    print(f'HTML report: {outfiles[i]}')
    print(f'JSON report: {json[i]}')
    
    #Run Fastplong
    #1. Trimming with fastplong
    subprocess.run(f'fastplong -i {infiles[i]} -o {reads[i]} -Q 0 >> {log} 2>&1', shell = True)
    #2. Move other fastplong outputs to the same folder
    subprocess.run(f'mv fastplong.html {outfiles[i]}', shell = True)
    subprocess.run(f'mv fastplong.json {json[i]}', shell = True)

