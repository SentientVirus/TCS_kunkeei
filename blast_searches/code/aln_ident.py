#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 17:52:07 2026

Script to calculate the percentage of identity and similar metrics between the 
RR-TF and reference sequences. This script is independent from the rest of the
pipeline, and hence isn't included in the Snakefile.'

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, re
import subprocess

# =============================================================================
# 1. Define regular expressions and functions
# =============================================================================
p = re.compile('\((.*)\)')
q = re.compile('Score: (.*)')

def run_mafft(infile, outfile, cores):
    command = f'mafft-linsi --thread {cores} {infile} > {outfile} 2> {log}'
    subprocess.run(command, shell = True)

# =============================================================================
# 2. Set paths to inputs and outputs
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/ugc00027/blast_searches' #Set working directory
log = f'{workdir}/logs/aln_indentity.log' #Set log file

for gene in ['RR-TF', 'HK']:
    infile = f'{workdir}/ref_aln/{gene}/{gene}.faa' #Set input file with all the RR-TF protein FASTA sequences
    out_mafft = infile.replace('.faa', '.mafft.faa') #Set the aligned output
    out_needle = infile.replace('.faa', '_needle.tsv') #Set the Needle output
    sep_dir = f'{workdir}/ref_aln/{gene}/separated_seqs' #Directory with the individual protein FASTA sequences
    in_prots = [f'{sep_dir}/{file}' for file in os.listdir(sep_dir) if file.endswith('.faa')] #Individual protein FASTA sequences
    t = 4 #No. of cores to use in MAFFT
        
    # =============================================================================
    # 3. Perform an alignment
    # =============================================================================
    
    run_mafft(infile, out_mafft, t)
    
    # =============================================================================
    # 4. Run EMBOSS Needle on the sequences and save results to files
    # =============================================================================
    
    with open(out_needle, 'w') as out: #Overwrite output file with global statistics
        out.write('seq1\tseq2\t%identity\t%similarity\t%gaps\tscore\n') #Write headers to file
        
    
    for i in range(len(in_prots)-1): #Loop through FASTA files with individual proteins
        aseq_name = os.path.basename(in_prots[i]).replace('.faa', '') #Retrieve the name of the first sequence
        for j in range(i + 1, len(in_prots)): #Loop again, excluding the first sequence and comparisons that were already performed
            bseq_name = os.path.basename(in_prots[j]).replace('.faa', '') #Retrieve the name of the second sequence
            # outfile = ... #This can replace stdout if an output file is wanted instead
            #Set the command to run EMBOSS Needle
            command = f'needle -aseq {in_prots[i]} -bseq {in_prots[j]} -outfile stdout -sprotein 1 -aformat simple -gapopen 10 -gapextend 0.5 2>> {log}'
            result = subprocess.run(command, shell = True, capture_output = True) #Run Needle and capture the output into a variable
            stdout_var = str(result.stdout).replace('\\n', '\n').split('\n') #Split the output by line breaks
            ident = p.search(stdout_var[25]).group(1).replace('%', '') #Retrieve the % identity
            simil = p.search(stdout_var[26]).group(1).replace('%', '') #Retrieve the % similarity
            gaps = p.search(stdout_var[27]).group(1).replace('%', '') #Retrieve the % gaps
            score = q.search(stdout_var[28]).group(1) #Retrieve the score
            with open(out_needle, 'a') as out: #Open the output file in append mode
                out.write(f'{aseq_name}\t{bseq_name}\t{ident}\t{simil}\t{gaps}\t{score}\n') #Write the pairwise Needle results to the file
            print(aseq_name, bseq_name) #Print the name of the proteins
            print(ident, simil, gaps, score) #Print the similarity metrics
    
    
        