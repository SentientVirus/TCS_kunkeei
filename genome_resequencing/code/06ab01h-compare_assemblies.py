# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:52:31 2026

Script to Blast the newly-sequenced genomes against each other.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
import logging, traceback
import sys
import subprocess

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0]

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger = logging.getLogger() #Create a logger
    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Defining inputs
# =============================================================================

infiles = sorted(snakemake.input) #Input assemblies
outfiles = sorted(snakemake.output) #Path to outputs

outpath = os.path.dirname(outfiles[0]) #Output directory

if not os.path.exists(outpath): #If the output directory does not exist
    os.makedirs(outpath) #Create it

# =============================================================================
# 2. Running Blast       
# =============================================================================

n = 0 #Counting variable to retrieve the right output file
for i in range(0, len(infiles)-1): #Loop through input files
    query = infiles[i] #Set query
    isolate1 = os.path.basename(query).replace('.fasta', '') #Set isolate name of query
    for j in range(i+1, len(infiles)): #Loop through input files, starting by the one after the query
        subject = infiles[j] #Set subject
        isolate2 = os.path.basename(subject).replace('.fasta', '') #Set isolate name of subject
        outfile = outfiles[n] #Set output file
        print(isolate1, isolate2, outfile) #Print isolate names and output file
        S = f'blastn -query {query} -out {outfile} -subject {subject} -outfmt 0' #Define BLAST command to run
        subprocess.run(S, shell = True) #Run the command
        n += 1 #Increase counting variable by one
