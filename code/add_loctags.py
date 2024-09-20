#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:41:24 2024

Script to add locus tags to the proteomics outputs

@author: Marina Mota Merlo
"""
import os
import pandas as pd
from Bio import GenBank

gbff = os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff'
infile = os.path.expanduser('~') + '/proteomics/files/MaxQuant_results.tsv'
# infile1 = os.path.expanduser('~') + '/proteomics/files/up_downregulated.tsv'
# infile2 = os.path.expanduser('~') + '/proteomics/files/unique_proteins.tsv'
# infiles = [infile1, infile2]


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
                

def add_loctags(df, rep_dict, out_file):
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
    
    outdir = os.path.dirname(out_file)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    df.to_csv(out_file, sep = '\t', index = False)
    
# for infile in infiles:
#     df = pd.read_csv(infile, sep = '\t')

#     df['Locus tags'] = df.loc[:, 'Protein IDs']
#     df['Majority locus tags'] = df.loc[:, 'Majority protein IDs']
    
#     outfile = infile.replace('files', 'files/loci')
    
#     add_loctags(df, replace_dict, outfile)

df = pd.read_csv(infile, sep = '\t')

df['Locus tags'] = df.loc[:, 'Protein IDs']
df['Majority locus tags'] = df.loc[:, 'Majority protein IDs']

outfile = infile.replace('files', 'files/loci')

add_loctags(df, replace_dict, outfile)