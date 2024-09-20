#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 13:33:52 2024

@author: marina
"""
import biolib
import os

infile = os.path.expanduser('~') + '/Akunkeei_files/faa/H3B1-04J_protein.faa'
outdir = os.path.expanduser('~') + '/proteomics/results/DeepTMHMM'

if not os.path.exists(outdir):
    os.makedirs(outdir)


deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24') # Load the version of DeepTMHMM that can run locally

biolib.utils.STREAM_STDOUT = True # Print progress
deeptmhmm_job = deeptmhmm.cli(args = f'--fasta {infile}', machine = 'local') # Run locally
deeptmhmm_job.save_files(outdir) # Set output direcory
