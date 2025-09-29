#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:11:46 2025

When the InterProScan run finishes, I should create another script that
adds the locus tags of the genes, instead of just the protein IDs. I can use
the one for Project 1 and modify it.

Once I have that, I should make another script to plot the predicted domains
as rectangles over the gene (not with Ete3, because there is no tree, but
with good old matplotlib)

This script retrieves the genes that are equivalent to the adhesins from
H3B1-04J from the A. kunkeei genomes.

@author: Marina Mota-Merlo
"""

import os, pandas as pd
import re

def is_sublist(sublist, main_list):
    sum_list = sum([element in main_list for element in sublist])
    if sum_list == len(sublist):
        return True
    else: return False

outfile = os.path.expanduser('~') + '/adhesins/results/adhesin_list.tsv'
annot_dir = os.path.expanduser('~') + '/adhesins/interproscan/locus_tags'

if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))

interpro_files = [file for file in os.listdir(annot_dir) if file.endswith('.tsv')]

pfam_domains = ['PF06458', 'PF19087', 'PF05737', #MucBP, DUF5776, collagen-binding
                'PF13632', 'PF17966', 'PF17965',  #gtf2, MucB2, MucBP_2
                'PF00746'] #LPXTG   'PF19258', 'TIGR03715' KxYKxGKxW SP

pfam_names = ['MucBP', 'DUF5776', 'collagen-binding', 'Gtf2', 'MucB2', 
              'MucBP_2', 'LPXTG']

genes_domains = {}

for file in interpro_files:
    with open(f'{annot_dir}/{file}') as ips:
        annot_df = pd.read_csv(ips, sep = '\t', header = None)
        pfam_df = annot_df[annot_df[4].isin(pfam_domains)]
        
        for index, row in pfam_df.iterrows():
            if '3M0' not in row[0] and row[0] not in genes_domains.keys():
                genes_domains[row[0]] = [[row[4]], [row[5]]]
            elif '3M0' not in row[0] and row[4] not in genes_domains[row[0]][0]:
                genes_domains[row[0]][0].append(row[4])
                genes_domains[row[0]][1].append(row[5])
        #Filter out the lines that are not Pfam
        #In the Pfam lines, search for any of the domains
        #If any of the domains is in the protein, record this information in a dictionary
        
tag_00510 = []
tag_00520 = []
tag_01020 = []
tag_adh = []
tag_01250 = []
tag_14310 = []
adh_dict = {}
strains = []
for key in list(genes_domains.keys()):
    strain = key.split('_')[0].replace('AKU', '').replace('AAP', '').replace('K2W83', 'DSMZ').replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('MUB42', 'HNS-8').replace('VQ058', 'GYUN-333')
    if re.search('H[0-9]', strain):
        strain = strain[:4] + '-' + strain[4:]
    if strain not in strains:
        strains.append(strain)
    if strain not in adh_dict.keys():
        adh_dict[strain] = ['', '']
    if genes_domains[key][0] == ['PF05737']: #Collagen-binding
        tag_00510.append(key)
    elif genes_domains[key][0] == ['PF13632']: #Gtf2
        tag_00520.append(key)
    elif is_sublist(sorted(['PF00746', 'PF17966', 'PF17965']), sorted(genes_domains[key][0])): #LPXTG + mucin-binding
        if 'PLPX' in key or int(key.replace('RS', '').split('_')[1]) < 1200 and key != 'K2W83_RS00655':
            tag_01020.append(key)
            adh_dict[strain][0] += f'{key}, '
        else:
            tag_adh.append(key)
            adh_dict[strain][1] += f'{key}, '
            
    elif genes_domains[key][0] == ['PF19087']: # and key.split('_')[0] not in ['K2W83', 'LDX55'] and int(key.split('_')[1]) < 1600:
        tag_01250.append(key)
    elif genes_domains[key][0] == ['PF19087']: #and key.split('_')[0] in ['K2W83', 'LDX55'] and int(key.split('_')[1].replace('RS', '')) < 800:
        tag_01250.append(key)
    elif sorted(genes_domains[key][0]) == sorted(['PF17966', 'PF00746']): #MubB2 + LPXTG
        tag_14310.append(key)

for strain in strains:
    adh_dict[strain][0] = adh_dict[strain][0][:-2]
    adh_dict[strain][1] = adh_dict[strain][1][:-2]
    
    if adh_dict[strain][0] == '':
        adh_dict[strain][0] = '-'
    if adh_dict[strain][1] == '':
        adh_dict[strain][1] = '-'

    if adh_dict[strain] == ['-', '-']:
        del adh_dict[strain]

with open(outfile, 'w') as handle:
    handle.write('Strain\tAdhesin_1\tAdhesin_2\n')
    for key in adh_dict.keys():
        handle.write(f'{key}\t{adh_dict[key][0]}\t{adh_dict[key][1]}\n')