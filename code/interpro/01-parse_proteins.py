#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:15:16 2026

Script to divide the file with all the proteins from strain H3B1-04J into
several files, each with a maximum of 100 proteins, to run the InterProScan
online server.

@author: Marina Mota-Merlo
"""

import os
from Bio import SeqIO

gbk = os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff'
faa = gbk.replace('genomic', 'protein').replace('gbff', 'faa')
workdir = os.path.expanduser('~') + '/mucoid_project/snpseq00064'
outdir = f'{workdir}/data/interpro'

if not os.path.exists(outdir):
    os.makedirs(outdir)

loctag_dict = {}
with open(gbk) as handle:
    for record in SeqIO.parse(handle, 'genbank'):
        record_dict = {feature.qualifiers['protein_id'][0]: feature.qualifiers['locus_tag'][0] for feature in record.features if ('locus_tag' in feature.qualifiers.keys()) and ('protein_id' in feature.qualifiers.keys())}
        loctag_dict = record_dict | loctag_dict
        
i = 1
n = 1
to_write = []

with open(faa) as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        record.description = record.description.replace(f'{record.id} ', '')
        record.description = record.description.replace(f'{loctag_dict[record.id]} ', '')
        record.id += f'/{loctag_dict[record.id]}'
        # print(record.id)
        if i > 100:
            n += 1
            i = 1
            
        print(n, i)
            
        outfile = f'{outdir}/H3B1-04J_proteins{n}.faa'
        if i == 1:
            with open(outfile, 'w') as faa_out:
                SeqIO.write(record, faa_out, 'fasta')
        else:
            with open(outfile, 'a') as faa_out:
                SeqIO.write(record, faa_out, 'fasta')
        i += 1