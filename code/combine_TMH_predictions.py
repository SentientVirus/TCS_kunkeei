#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 16:22:43 2024

Script to compare the results of Phobius to those of DeepTMHMM and keep
the overlapping transmembrane helices, then save them to a file together
with the protein ID and protein locus tag.

@author: Marina Mota Merlo
"""
import os
import pandas as pd
from Bio import GenBank

gbff = os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff'
phob = os.path.expanduser('~') + '/proteomics/results/Phobius/H3B1-04J_phobius.txt'
DTMH = os.path.expanduser('~') + '/proteomics/results/DeepTMHMM/TMRs.gff3'
pred_out = os.path.expanduser('~') + '/proteomics/results/TMH_predictions/TMH.tab'

infile = os.path.expanduser('~') + '/proteomics/files/loci/MaxQuant_results.tsv'
outfile = os.path.expanduser('~') + '/proteomics/files/loci/MaxQuant_results_TMH.tsv'

if not os.path.exists(os.path.dirname(pred_out)):
    os.makedirs(os.path.dirname(pred_out))
    
with open(outfile, 'w') as out_pred:
    out_pred.write('locus_tag\tprotein_id\tstart\tend\n')

class TM:
    def __init__(self, prot_id, locus_tag, start, end):
        self.prot_id = prot_id
        self.locus_tag = locus_tag
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start
    def __str__(self):
        return f'TMH in protein {self.prot_id} ({self.start}-{self.end})'
    def write_TM(self, out_file):
        with open(out_file, 'a') as handle: 
            handle.write(f'{self.locus_tag}\t{self.prot_id}\t{self.start}\t{self.end}\n')

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

phobius_TM = []
with open(phob) as phobius:
    print('Loading Phobius results...')
    for line in phobius:
        if line.startswith('ID'):
           line = line.strip()
           prot_id = line.split(' ')[-1]
           loctag = replace_dict[prot_id]
           
        elif 'TRANSMEM' in line:
            line = line.strip()
            pos = line.split(' ')
            while '' in pos:
                pos.remove('')
            helix = TM(prot_id, loctag, pos[2], pos[3])
            print(helix)
            phobius_TM.append(helix)
                
deep_TM = []
with open(DTMH) as deepTM:
    print('Loading DeepTMHMM results...')
    for line in deepTM:
        if line.startswith('CAI') and 'TMhelix' in line:
            line = line.strip()
            line = line.split('\t')
            while '' in line:
                line.remove('')
                
            helix = TM(line[0], replace_dict[line[0]], line[2], line[3])
            print(helix)
            deep_TM.append(helix)
 
add_info = {}
for TMH1 in phobius_TM:
    for TMH2 in deep_TM:
        if TMH1.prot_id == TMH2.prot_id and (TMH1.start < TMH2.end and TMH1.end > TMH2.start):
            new_helix = TM(TMH1.prot_id, TMH1.locus_tag, max(TMH1.start, TMH2.start), min(TMH1.end, TMH2.end))
            new_helix.write_TM(pred_out)
            if new_helix.prot_id not in add_info.keys():
                add_info[new_helix.prot_id] = f'{new_helix.start}-{new_helix.end}'
            else:
                add_info[new_helix.prot_id] += f';{new_helix.start}-{new_helix.end}'
            

df = pd.read_csv(infile, sep = '\t')

df['TMH?'] = df.loc[:, 'Protein IDs']
df['TMH positions'] = df.loc[:, 'Protein IDs']
for index, row in df.iterrows():
    if df.loc[index, 'Majority protein IDs'] in add_info.keys():
        df.loc[index, 'TMH?'] = 'Y'
        df.loc[index, 'TMH positions'] = add_info[df.loc[index, 'Majority protein IDs']]
    else:
        df.loc[index, 'TMH?'] = 'N'
        df.loc[index, 'TMH positions'] = float('nan')
        
        
cols = list(df)

cols.insert(7, cols.pop(cols.index('TMH?')))
cols.insert(8, cols.pop(cols.index('TMH positions')))
df = df.loc[:, cols]

df.to_csv(outfile, sep = '\t', index = False)
        



