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

log = snakemake.log[0]

with open(log, 'w') as logfile: #Overwrite log file
    logfile.write('')
            
#Redirect stdout and stderr to log file
sys.stdout = open(log, 'a')
sys.stderr = open(log, 'a')

#Format the logging
logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================

infile = snakemake.input #The protein.faa file for strain H3B1-04J
outdir = snakemake.params.outdir #Path to outputs

logging.info(f'Create output directory {outdir} if it doesn\'t exist...')
if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

# =============================================================================
# 2. Run DeepTMHMM
# =============================================================================

logging.info('Load DeepTMHMM...')
deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24') #Load the version of DeepTMHMM that can run locally

logging.info(f'Run DeepTMHMM on {infile}...')
biolib.utils.STREAM_STDOUT = True #Print progress
deeptmhmm_job = deeptmhmm.cli(args = f'--fasta {infile}', machine = 'local') #Run DeepTMHMM locally on the input file
deeptmhmm_job.save_files(outdir) #Set output direcory
logging.info('Done!')