# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:52:31 2026

Script to Blast the newly-sequenced genomes against each other.

@author: Marina Mota-Merlo
"""
import os
import logging, traceback
import sys
import subprocess

# =============================================================================
# Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/ugc00027'
log = f'{workdir}/logs/test_blastn.log' #snakemake.log[0]

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
# Defining inputs
# =============================================================================

isolates = ['01', '02', '09', '10']

infiles = [f'{workdir}/assemblies/combined/{isol}.fasta' for isol in isolates]
outpath = f'{workdir}/test/blastn' #snakemake.params.outpath #Path to outputs
threads = 12 #snakemake.threads #Number of threads to be used

if not os.path.exists(outpath):
    os.makedirs(outpath)

# =============================================================================
# Running Blast       
# =============================================================================

for i in range(0, len(infiles)-1): #Loop through input files
    query = infiles[i]
    isolate1 = os.path.basename(query).replace('.fasta', '')
    for j in range(i+1, len(infiles)):
        subject = infiles[j]
        isolate2 = os.path.basename(subject).replace('.fasta', '')
        outfile = f'{outpath}/{isolate1}_vs_{isolate2}.txt'
        print(outfile)
        S = f'blastn -query {query} -out {outfile} -subject {subject} -outfmt 0' #Define BLAST command to run
        subprocess.run(S, shell = True) #Run the command
    
