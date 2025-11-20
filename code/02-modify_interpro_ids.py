#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 12:35:32 2025.
Modified from a script created by Marina Mota-Merlo on Thu Apr 13 11:45:45 2023

Script to change the column with protein IDs in the Interproscan output to 
a column with locus tags, so that these files can be navigated using the same
scripts as for the A. kunkeei annotations made by Dyrhage et al. (2022).

@author: Marina Mota-Merlo
"""

# =============================================================================
# Script to change protein IDs in Interproscan annotations to locus tags
# =============================================================================

import os
import pandas as pd
from Bio import SeqIO
import multiprocessing
from functools import partial
import time

start_time = time.time() 
count = 0

workdir = os.path.expanduser('~') + '/mucoid_project/adhesins/interproscan'
outpath = f'{workdir}/locus_tags'
path = os.path.expanduser('~') + '/Akunkeei_files/gbff'
new_strain_path = path.replace('gbff', 'new_genomes/gbff')
prot_loctag_out = os.path.expanduser('~') + '/mucoid_project/adhesins/metadata/prot_id_loctag.tsv'
infiles = [f'{path}/{file}' for file in os.listdir(path) if not os.path.isdir(f'{path}/{file}')] + [f'{new_strain_path}/{file}' for file in os.listdir(new_strain_path)]
no_cores = 36
prot2loctag = {}

if not os.path.exists(outpath):
    os.makedirs(outpath)
    
if not os.path.exists(os.path.dirname(prot_loctag_out)):
    os.makedirs(os.path.dirname(prot_loctag_out))

def process_infiles(file):
    prot2loctag = {}
    with open(file) as genb:
        info = SeqIO.parse(genb, 'genbank')
        strain = os.path.basename(file).split('_')[0]
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
        df2 = df.replace({0: prot2loctag})
        df2.to_csv(f'{outpath}/{strain}.tsv', sep = '\t', header = False, index = False)
    return prot2loctag

#Multithreading
if __name__ == '__main__': 
    pool = multiprocessing.Pool() 
    pool = multiprocessing.Pool(processes=no_cores)
    prot2loctag = {}
    for file in infiles:
        count += 1
        print(f'Processing strain {count}/{len(infiles)} ({os.path.basename(file).split("_")[0]})...')
        outputs = pool.map(partial(process_infiles), infiles)

with open(prot_loctag_out, 'w') as handle:
    handle.write('Locus_tag\tProtein_id\n')
    [{handle.write(f'{value}\t{key}\n') for (key, value) in output.items()} for output in outputs]
    
end_time = time.time() - start_time
print(f'This script took {end_time:2f} with {no_cores} processes')
