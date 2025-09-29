#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 12:35:32 2025.
Modified from a script created on Thu Apr 13 11:45:45 2023

Script to change the column with protein IDs in the Interproscan output to 
a column with locus tags, so that these files can be navigated using the same
scripts as for the A. kunkeei annotations made by Dyrhage et al. (2022).

@author: Marina Mota Merlo
"""

# =============================================================================
# Script to change protein IDs in Interproscan annotations to locus tags
# =============================================================================

import os
import pandas as pd
from Bio import SeqIO

workdir = os.path.expanduser('~') + '/adhesins/interproscan'
outpath = os.path.expanduser('~') + '/adhesins/interproscan/locus_tags'
path = os.path.expanduser('~') + '/Akunkeei_files/gbff'
new_strain_path = path.replace('gbff', 'new_genomes/gbff')
infiles = [f'{path}/{file}' for file in os.listdir(path) if not os.path.isdir(f'{path}/{file}')] + [f'{new_strain_path}/{file}' for file in os.listdir(new_strain_path)]

prot2loctag = {}

if not os.path.exists(outpath):
    os.makedirs(outpath)

for file in infiles:
# for s in new_strains:
    strain = os.path.basename(file).replace('_genomic.gbff', '')
    with open(file) as genb:
        info = SeqIO.parse(genb, 'genbank')
        for record in info:
            for feature in record.features:
                tag = feature.qualifiers.get('locus_tag')
                prot_id = feature.qualifiers.get('protein_id')
                if tag and prot_id:
                    tag = tag[0]
                    prot_id = prot_id[0]
                    prot2loctag[prot_id] = tag
    with open(f'{workdir}/{strain}.tsv') as df_file:
        df = pd.read_csv(df_file, sep = '\t', header = None)
        protein_ids = df[0]
        df2 = df.replace({0: prot2loctag})
        df2.to_csv(f'{outpath}/{strain}.tsv', sep = '\t', header = False, index = False)
