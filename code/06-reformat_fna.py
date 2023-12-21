#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:52:44 2023

@author: marina

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

input_files = snakemake.input
output_files = snakemake.output

line_length = 80


# =============================================================================
# Changing length of lines in fasta files
# =============================================================================

for i in range(len(input_files)):
    input_file = input_files[i]
    output_file = output_files[i]
    print(input_file, output_file)
    with open(f'{input_file}') as f_in, open(f'{output_file}', 'w') as f_out:
        sequence = ''
        n = 1
        for line in f_in:
            if line.startswith('>'):
                line = f'>chromosome_{n}\n'
                n += 1
                if sequence:
                    # Write the sequence with increased line length
                    f_out.write('\n'.join([sequence[i:i+line_length] for i in range(0, len(sequence), line_length)]) + '\n')
                    sequence = ''
                f_out.write(line)
            else:
                sequence += line.strip()
        
        # Write the last sequence with increased line length
        if sequence:
            f_out.write('\n'.join([sequence[i:i+line_length] for i in range(0, len(sequence), line_length)]))
            f_out.write('\n')
