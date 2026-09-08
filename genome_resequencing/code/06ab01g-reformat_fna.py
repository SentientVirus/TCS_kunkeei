#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:52:44 2023

Script to change the format of nucleotide FASTA files so that they show 80
nucleotides per line instead of 60.

@author: Marina Mota-Merlo

"""
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
# Defining inputs
# =============================================================================

input_files = snakemake.input #Path to input files
output_files = snakemake.output #Path to output files

line_length = 80 #Desired length of lines in the FASTA file


# =============================================================================
# Changing length of lines in fasta files
# =============================================================================

for i in range(len(input_files)): #Loop through input files
    input_file = input_files[i] #Get path to ith input file
    output_file = output_files[i] #Get path to ith output file
    print(input_file, output_file) #Print both
    with open(f'{input_file}') as f_in, open(f'{output_file}', 'w') as f_out: #Open input and output files
        sequence = '' #Create string to store the sequence
        n = 1 #Initialize count variable
        for line in f_in: #Loop through lines in the input file
            if line.startswith('>'): #If the line is a FASTA header
                line = f'>chromosome_{n}\n' #Change the text in the header
                n += 1 #Increase count variable
                if sequence: #If the sequence string is not empty
                    # Write the sequence with increased line length
                    f_out.write('\n'.join([sequence[i:i+line_length] for i in range(0, len(sequence), line_length)]) + '\n') #Write the sequence with 80 nucleotides per line
                    sequence = '' #Re-initialize the sequence variable
                f_out.write(line) #Write the header line to output
            else: #If the line is not a FASTA header
                sequence += line.strip() #Store it in the string
        
        # Write the last sequence with increased line length
        if sequence:
            f_out.write('\n'.join([sequence[i:i+line_length] for i in range(0, len(sequence), line_length)]))
            f_out.write('\n')
