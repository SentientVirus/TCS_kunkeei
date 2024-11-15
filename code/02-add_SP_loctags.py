#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:41:24 2024

Script to add locus tags and SignalP predictions to the proteomics outputs

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================
import os
import pandas as pd
from Bio import GenBank

# =============================================================================
# 0. Define functions
# =============================================================================

def add_loctags(df, rep_dict):
    df['Locus tags'] = df.loc[:, 'Protein IDs']
    df['Majority locus tags'] = df.loc[:, 'Majority protein IDs']
    for index, row in df.iterrows():
        id_list = df.loc[index, 'Locus tags'].split(';')
        locus_list = []
        
        for protid in id_list:
            for prot in replace_dict.keys():
                rep = df.loc[index, 'Majority locus tags']
                if prot in rep:
                    df.loc[index, 'Majority locus tags'] = rep.replace(prot, replace_dict[prot])
                if prot in protid:
                    locus = protid.replace(prot, replace_dict[prot])
                    
            
            locus_list.append(locus)
        
        locus_tags = ';'.join(locus_list)
        
        df.loc[index, 'Locus tags'] = locus_tags
        
        
        
    cols = list(df)
    
    cols.insert(0, cols.pop(cols.index('Locus tags')))
    cols.insert(1, cols.pop(cols.index('Majority locus tags')))
    
    df = df.loc[:, cols]
    
    return df

# =============================================================================
# 1. Define paths to inputs and outputs
# =============================================================================
workdir = os.path.expanduser('~') + '/proteomics'
gbff = os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff'
indir = f'{workdir}/files/parsed'
signalP = os.path.expanduser('~') + '/proteomics/results/SignalP/H3B1-04J_SignalP.txt'
SP_tab = signalP.replace('txt', 'tsv')
outdir = f'{workdir}/files/loci'

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
# =============================================================================
# 2. Add information to the input file and write it to an output
# =============================================================================
infiles = [file for file in os.listdir(indir) if file.endswith('.tsv')]
for infile in infiles:
    outfile = f'{outdir}/{infile}'
    
    replace_dict = {}
    with open(gbff) as handle:
        for record in GenBank.parse(handle):
            for feature in record.features:
                loctag = ''
                prot_id = ''
                for qual in feature.qualifiers:
                    if 'locus_tag' in qual.key:
                        loctag = qual.value.strip('"')
                    elif 'protein_id' in qual.key:
                        prot_id = qual.value.strip('"')
                if loctag != '' and prot_id != '':
                    replace_dict[prot_id] = loctag
                    print(f'Saving protein {prot_id} with locus tag {loctag}...')
    
    df = pd.read_csv(f'{indir}/{infile}', sep = '\t')
    
    df = add_loctags(df, replace_dict)
    
    with open(signalP) as SP_file, open(SP_tab, 'w') as tabfile:
        for line in SP_file:
            line = line.split(' ')
            while '' in line:
                line.remove('')
            if 'name' in line:
                line = line[1:]
                line[4] = 'end'
            line_text = '\t'.join(line) + '\n'
            if '#' not in line_text:
                tabfile.write(line_text)
            
            
    SP_df = pd.read_csv(SP_tab, sep = '\t')
    
    df['SP?'] = df.loc[:, 'Protein IDs']
    df['SP positions'] = df.loc[:, 'Protein IDs']
    for index, row in df.iterrows():
        for index2, row2 in SP_df.iterrows():
            if df.loc[index, 'Majority protein IDs'] == SP_df.loc[index2, 'name']:
                print(f'Adding signal peptide prediction to protein {df.loc[index, "Majority protein IDs"]}')
                df.loc[index, 'SP?'] = SP_df.loc[index2, '?']
                if SP_df.loc[index2, '?'] == 'Y':
                    df.loc[index, 'SP positions'] = f'1-{int(SP_df.loc[index2, "end"])-1}'
                else:
                    df.loc[index, 'SP positions'] = float('nan')
                
    cols = list(df)
    
    cols.insert(5, cols.pop(cols.index('SP?')))
    cols.insert(6, cols.pop(cols.index('SP positions')))
    df = df.loc[:, cols]
    
    df.to_csv(outfile, sep = '\t', index = False)
