#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 13:33:52 2024

Script to run DeepTMHMM.

OBS! This script prompts an error message when run with Snakemake and the
--use-conda option. Please, install and activate the python_env.yml
envrionment before running the script.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries
# =============================================================================

import biolib
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
# 1. Define inputs and outputs
# =============================================================================

infile = snakemake.input #os.path.expanduser('~') + '/Akunkeei_files/faa/H3B1-04J_protein.faa'
outdir = snakemake.params.outdir #os.path.expanduser('~') + '/proteomics/results/DeepTMHMM'

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

# =============================================================================
# 2. Run DeepTMHMM
# =============================================================================

deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24') #Load the version of DeepTMHMM that can run locally

biolib.utils.STREAM_STDOUT = True #Print progress
deeptmhmm_job = deeptmhmm.cli(args = f'--fasta {infile}', machine = 'local') #Run DeepTMHMM locally on the input file
deeptmhmm_job.save_files(outdir) #Set output direcory
