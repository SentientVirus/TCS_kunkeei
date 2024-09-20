#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 13:33:52 2024

@author: marina
"""
import biolib
import os

infile = os.path.expanduser('~') + '/Akunkeei_files/faa/H3B1-04J_protein.faa'

deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24')

biolib.utils.STREAM_STDOUT = True # Stream progress from app in real time
deeptmhmm_job = deeptmhmm.cli(args = f'--fasta {infile}', machine = 'local') # Blocks until done
deeptmhmm_job.save_files('results') # Saves all results to `result` dir
