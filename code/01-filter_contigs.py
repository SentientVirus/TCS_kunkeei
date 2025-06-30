#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:26:50 2023

Script to reverse the strand of the genome of isolate 10, which has the
reverse strand as forward strand in the raw assembly.

@author: Marina Mota-Merlo
"""

import os
from Bio import SeqIO
import logging, traceback
from Bio.SeqUtils import GC
import pandas as pd

# =============================================================================
# Logging
# =============================================================================

#logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
#                    format = '%(asctime)s %(message)s',
#                    datefmt = '%Y-%m-%d %H:%M:%S')

#def handle_exception(exc_type, exc_value, exc_traceback):
#    if issubclass(exc_type, KeyboardInterrupt):
#        sys.__excepthook__(exc_type, exc_value, exc_traceback)
#        return

#    logger.error(''.join(["Uncaught exception: ",
#                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
#                          ]))

#sys.excepthook = handle_exception

#sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# Define inputs and outputs
# =============================================================================

indir = '../assemblies'
infiles = [f'{file[0]}/{file[2][2]}' for file in os.walk(indir) if file[0].startswith(f'{indir}/') and len(file[0]) == len(indir) + 3] # if file.endswith('assembly.fasta')
assemblies =  [f'{file[0]}/{file[2][1]}' for file in os.walk(indir) if file[0].startswith(f'{indir}/') and len(file[0]) == len(indir) + 3]
outfiles = [f'{file[0]}/{file[2][1].replace(".fasta", "_filtered.fasta")}' for file in os.walk(indir) if file[0].startswith(f'{indir}/') and len(file[0]) == len(indir) + 3]

print(infiles)

# =============================================================================
# Reverse strand and save to file
# =============================================================================

for i in range(len(infiles)):
    to_exclude = []
    with open(infiles[i]) as handle:
        print('Filename: ', infiles[i])
        df = pd.read_csv(handle, sep = '\t')
        for index, contig in df.iterrows():
            contig_name = contig['#seq_name']
            print('Contig name: ', contig_name)
            if contig['cov.'] < 30:
                to_exclude.append(contig_name)
        print('Contigs to exclude: ', to_exclude)
        
    with open(assemblies[i]) as assembly:
        records = SeqIO.parse(assembly, 'fasta')
        to_write = []
        for record in records:
            if record.id not in to_exclude:
                to_write.append(record)
                
    with open(outfiles[i], 'w') as outfile:
        SeqIO.write(to_write, outfile, 'fasta')
                
    # records = SeqIO.parse(file, 'fasta')
    # for record in records:
    #     print('Filename: ', file)
    #     print('Contig name: ', record.id)
    #     print('GC content: ', GC(record.seq))
    #     print('Contig length: ', len(record.seq))
    #     AT_ratio = str(record.seq).count('A')/str(record.seq).count('T')
    #     GC_ratio = str(record.seq).count('G')/str(record.seq).count('C')
    #     print('Nucleotide proportions (A/T & G/C): ', f'{AT_ratio:2f} ', f'{GC_ratio:2f}')
    #     if not (0.75 <= AT_ratio <= 1/0.75 and 0.75 <= GC_ratio <= 1/0.75):
    #         to_exclude.append(record.id)
    # print('Contigs to exclude: ', to_exclude)

#with open(filename) as reverse_stranded:
#    records = SeqIO.parse(filename, 'fasta')
#    for record in records:
#        record.seq = record.seq.reverse_complement()
#        record_list.append(record)

#with open(outfile, 'w') as output:
#    SeqIO.write(record_list, output, 'fasta')
