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
from Bio import SeqIO

def is_sublist(sublist, main_list):
    sum_list = sum([element in main_list for element in sublist])
    if sum_list == len(sublist):
        return True
    else: return False

outfile = os.path.expanduser('~') + '/adhesins/results/adhesin_list.tsv'
annot_dir = os.path.expanduser('~') + '/adhesins/interproscan/locus_tags'
faa_dir = os.path.expanduser('~') + '/Akunkeei_files/faa'
faa_outfile = os.path.expanduser('~') + '/adhesins/sequences/Muc_adhesins.faa'
faa_repset = os.path.expanduser('~') + '/adhesins/sequences/Muc_adhesins_repset.faa'
prot2loctag = os.path.expanduser('~') + '/adhesins/metadata/prot_id_loctag.tsv'

repr_strains = ['DSMZ12361', 'IBH001', 'GYUN-333', 'HNS-8', 'A0901', 
                'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805', 
                'Fhon2', 'G0102', 'G0403', 'H1B1-04J', 'H1B1-05A', 
                'H1B3-02M', 'H3B1-11M', 'H3B1-04J', 'H3B1-04X', 'H3B1-03M', #The last one is an extra strain, added because of the plasmid gene
                'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 'H4B2-06J', 'H4B4-02J', 
                'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 
                'H4B5-04J', 'H4B5-05J', 'MP2', 'Fhon13']

[os.makedirs(os.path.dirname(out)) for out in [outfile, faa_outfile] if not os.path.exists(os.path.dirname(out))]
# if not os.path.exists(os.path.dirname(outfile)):
#     os.makedirs(os.path.dirname(outfile))

interpro_files = [file for file in os.listdir(annot_dir) if file.endswith('.tsv')]

pfam_domains = ['PF06458', 'PF19087', 'PF05737', #MucBP, DUF5776, collagen-binding
                'PF13632', 'PF17966', 'PF17965',  #gtf2, MucB2, MucBP_2
                'PF00746'] #LPXTG   'PF19258', 'TIGR03715' KxYKxGKxW SP

pfam_names = ['MucBP', 'DUF5776', 'collagen-binding', 'Gtf2', 'MucB2', 
              'MucBP_2', 'LPXTG']

with open(prot2loctag) as handle:
    loctag_dict = {line.split('\t')[1].strip(): line.split('\t')[0] for line in handle}

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
to_retrieve = {}
for key in list(genes_domains.keys()):
    strain = key.split('_')[0].replace('AKU', '').replace('AAP', '').replace('K2W83', 'DSMZ').replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('MUB42', 'HNS-8').replace('VQ058', 'GYUN-333')
    if re.search('H[0-9]', strain):
        strain = strain[:4] + '-' + strain[4:]
    strain = strain.replace('FHON', 'Fhon').replace('DSMZ', 'DSMZ12361')
    if strain not in strains:
        strains.append(strain)
        to_retrieve[strain] = []
    if strain not in adh_dict.keys():
        adh_dict[strain] = ['', '']
    if genes_domains[key][0] == ['PF05737']: #Collagen-binding
        tag_00510.append(key)
    elif genes_domains[key][0] == ['PF13632']: #Gtf2
        tag_00520.append(key)
    elif is_sublist(sorted(['PF17966', 'PF17965']), sorted(genes_domains[key][0])) or key == 'AKUH3B104X_PLPX00300': #LPXTG + mucin-binding
        if ('PLPX' in key or int(key.replace('RS', '').split('_')[1]) < 1200 or 'MUB' in key) and key not in ['K2W83_RS00655', 'AAPFHON13_01040']:
            tag_01020.append(key)
            adh_dict[strain][0] += f'{key}, '
        else:
            tag_adh.append(key)
            adh_dict[strain][1] += f'{key}, '
        to_retrieve[strain].append(key)
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
    
    if to_retrieve[strain] == []:
        del to_retrieve[strain]
        

with open(outfile, 'w') as handle:
    handle.write('Strain\tAdhesin_locus1\tAdhesin_locus2\n')
    for key in adh_dict.keys():
        handle.write(f'{key}\t{adh_dict[key][0]}\t{adh_dict[key][1]}\n')
        
with open(faa_outfile, 'w') as faa_out, open(faa_repset, 'w') as faa_rep:
    for genome in to_retrieve.keys():
        faa = f'{faa_dir}/{genome}_protein.faa'
        check = False
        if genome == 'HNS-8':
            faa = faa.replace('/faa', '/new_genomes/faa')
        if genome in repr_strains:
            check = True
        for record in SeqIO.parse(faa, 'fasta'):
            if loctag_dict[record.id] in to_retrieve[genome]:
                record.id = loctag_dict[record.id]
                record.description = ''
                SeqIO.write(record, faa_out, 'fasta')
                if check:
                    SeqIO.write(record, faa_rep, 'fasta')
            
        