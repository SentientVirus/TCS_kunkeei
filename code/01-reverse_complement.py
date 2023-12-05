#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:26:50 2023

@author: marina
"""

from Bio import SeqIO

filename = 'data/004.fna'
record_list = []

with open(filename) as reverse_stranded:
    records = SeqIO.parse(filename, 'fasta')
    for record in records:
        record.seq = record.seq.reverse_complement()
        record_list.append(record)

outfile = 'data/rev004.fna'        
with open(outfile, 'w') as output:
    SeqIO.write(record_list, output, 'fasta')
