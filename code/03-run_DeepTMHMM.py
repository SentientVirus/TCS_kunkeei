#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 13:33:52 2024

Script to run DeepTMHMM.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries
# =============================================================================

import biolib
import os

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================

infile = os.path.expanduser('~') + '/Akunkeei_files/faa/H3B1-04J_protein.faa'
outdir = os.path.expanduser('~') + '/proteomics/results/DeepTMHMM'

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

# =============================================================================
# 2. Run DeepTMHMM
# =============================================================================

deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24') #Load the version of DeepTMHMM that can run locally

biolib.utils.STREAM_STDOUT = True #Print progress
deeptmhmm_job = deeptmhmm.cli(args = f'--fasta {infile}', machine = 'local') #Run DeepTMHMM locally on the input file
deeptmhmm_job.save_files(outdir) #Set output direcory
