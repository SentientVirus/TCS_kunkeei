#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:11:46 2025

This script retrieves the genes that are equivalent to the adhesins from
H3B1-04J from the A. kunkeei genomes and creates a file with the locus tags
and another file with the sequences for each gene.

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

class adhesin:
    '''
    Class to store the locus tag, type and strain of each adhesin.
    '''
    
    def __init__(self, locus_tag, strain, adh_type = '?', ref_locus = '', pfams = [], domain_names = []): #Define class
        self.locus_tag = locus_tag
        self.ref_locus = ref_locus
        self.strain = strain
        self.type = adh_type
        self.pfams = pfams
        self.domain_names = domain_names
        
    def get_type(self):
        if self.pfams == ['PF05737']:
            self.type = 'collagen-binding'
            self.ref_locus = 'AKUH3B104J_00510'
        elif self.pfams == ['PF13632']:
            self.type = 'Gtf2'
            self.ref_locus = 'AKUH3B104J_00520'
        elif self.pfams == ['PF19087']:
            self.type = 'SH3b'
            self.ref_locus = 'AKUH3B104J_01250'
        elif sorted(self.pfams) == sorted(['PF17966', 'PF00746']):
            self.type = 'MubB2+LPXTG'
            self.ref_locus = 'AKUH3B104J_14310'
        elif is_sublist(sorted(['PF17966', 'PF17965']), self.pfams) or self.locus_tag == 'AKUH3B104X_PLPX00300':
            self.type = 'MucBP+LPXTG'
            if ('PLPX' in self.locus_tag or int(self.locus_tag.replace('RS', '').split('_')[1]) < 1200 or 'MUB' in self.locus_tag) and self.locus_tag not in ['K2W83_RS00655', 'AAPFHON13_01040']:
                self.ref_locus = 'AKUH3B104J_01020'
            else:
                self.ref_locus = 'None'
    def __str__(self): #What print(class) returns
        return f'<Adhesin {self.locus_tag} from strain {self.strain} is {self.type}>'


adhesins = ['MucBP+LPXTG', 'MubB2+LPXTG', 'Gtf2', 'collagen-binding', 'SH3b']
    
workdir = os.path.expanduser('~') + '/adhesins'
annot_dir = f'{workdir}/interproscan/locus_tags'
faa_dir = os.path.expanduser('~') + '/Akunkeei_files/faa'
prot2loctag = f'{workdir}/metadata/prot_id_loctag.tsv'
outdir = f'{workdir}/results'
seqdir = f'{workdir}/sequences'

pfam_domains = ['PF06458', 'PF19087', 'PF05737', #MucBP, DUF5776, collagen-binding
                'PF13632', 'PF17966', 'PF17965',  #gtf2, MucB2, MucBP_2
                'PF00746'] #LPXTG   'PF19258', 'TIGR03715' KxYKxGKxW SP

pfam_names = ['MucBP', 'DUF5776', 'collagen-binding', 'Gtf2', 'MucB2', 
              'MucBP_2', 'LPXTG']

repr_strains = ['DSMZ12361', 'IBH001', 'GYUN-333', 'HNS-8', 'A0901', 
                'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805', 
                'Fhon2', 'G0102', 'G0403', 'H1B1-04J', 'H1B1-05A', 
                'H1B3-02M', 'H3B1-11M', 'H3B1-04J', 'H3B1-04X', 'H3B1-03M', #The last one is an extra strain, added because of the plasmid gene
                'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 'H4B2-06J', 'H4B4-02J', 
                'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 
                'H4B5-04J', 'H4B5-05J', 'MP2', 'Fhon13']

interpro_files = [f'{annot_dir}/{file}' for file in os.listdir(annot_dir) if file.endswith('.tsv') and 'M-0' not in file]

[os.makedirs(out) for out in [outdir, seqdir] if not os.path.exists(out)]

with open(prot2loctag) as handle:
    loctag_dict = {line.split('\t')[1].strip(): line.split('\t')[0] for line in handle}
    
# genes_domains = {}
gene_list = []
curr_adh = adhesin('', '')
adh_dict = {}
strains = []
to_retrieve = {}
for file in interpro_files:
    with open(file) as ips:
        annot_df = pd.read_csv(ips, sep = '\t', header = None)
        pfam_df = annot_df[annot_df[4].isin(pfam_domains)] #Filter out the lines that are not 
        strain = os.path.basename(file).split('_')[0].replace('.tsv', '')
        strains.append(strain)
        adh_dict[strain] = []
        
        for index, row in pfam_df.iterrows(): #In the Pfam lines, search for any of the domains
            if  row[0] not in gene_list:
                #If any of the domains is in the protein, record this information in an object of class adhesin
                # genes_domains[row[0]] = [[row[4]], [row[5]]]
                gene_list.append(row[0])
                next_adh = adhesin(row[0], strain, pfams = [])
                print(next_adh)
                if next_adh.locus_tag != curr_adh.locus_tag:
                    curr_adh.get_type()
                    if curr_adh.locus_tag != '':
                        adh_dict[curr_adh.strain].append(curr_adh)
                    curr_adh = next_adh
            curr_adh.pfams.append(row[4])
            curr_adh.domain_names.append(row[5])
                
for adhesin in adhesins:
    outfile = f'{outdir}/{adhesin}_list.tsv'
    faa_outfile = f'{seqdir}/{adhesin}.faa'
    faa_repset = f'{seqdir}/{adhesin}_repset.faa'
    with open(outfile, 'w') as handle:
        if adhesin == 'MucBP+LPXTG':
            handle.write('Strain\tAdhesin_locus1\tAdhesin_locus2\n')
        else:
            handle.write('Strain\tAdhesin_locus\n')
        for strain in sorted(strains):
            adh_type = [adh for adh in adh_dict[strain] if adh.type == adhesin]
            to_retrieve[strain] = [adh.locus_tag for adh in adh_type]
            if adhesin == 'MucBP+LPXTG':
                adh1 = [adh for adh in adh_type if adh.ref_locus == 'AKUH3B104J_01020']
                adh2 = [adh for adh in adh_type if adh.ref_locus == 'None']
                if len(adh1) == 0:
                    adh1_string = '-'
                else:
                    adh1_string = ''
                    for adh in adh1:
                        adh1_string += adh.locus_tag + ', '
                    adh1_string = adh1_string[:-2]
                if len(adh2) == 0:
                    adh2_string = '-'
                else:
                    adh2_string = ''
                    for adh in adh2:
                        adh2_string += adh.locus_tag + ', '
                    adh2_string = adh2_string[:-2]
                    
                if not (adh1_string == '-' and adh2_string == '-'): 
                    handle.write(f'{strain}\t{adh1_string}\t{adh2_string}\n')
            else:
                [handle.write(f'{strain}\t{adh.locus_tag}\n') for adh in adh_type]
                
    with open(faa_outfile, 'w') as faa_out, open(faa_repset, 'w') as faa_rep:
        for genome in sorted(list(to_retrieve.keys())):
            faa = f'{faa_dir}/{genome}_protein.faa'
            check = False
            if genome == 'HNS-8' or genome == 'GYUN-333':
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
                
            
                
                
        
        
# tag_00510 = []
# tag_00520 = []
# tag_01020 = []
# tag_adh = []
# tag_01250 = []
# tag_14310 = []
# adh_dict = {}
# strains = []
# to_retrieve = {}
# for key in list(genes_domains.keys()):
#     strain = key.split('_')[0].replace('AKU', '').replace('AAP', '').replace('K2W83', 'DSMZ').replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('MUB42', 'HNS-8').replace('VQ058', 'GYUN-333')
#     if re.search('H[0-9]', strain):
#         strain = strain[:4] + '-' + strain[4:]
#     strain = strain.replace('FHON', 'Fhon').replace('DSMZ', 'DSMZ12361')
#     if strain not in strains:
#         strains.append(strain)
#         to_retrieve[strain] = []
#     if strain not in adh_dict.keys():
#         adh_dict[strain] = ['', '']
#     if genes_domains[key][0] == ['PF05737']: #Collagen-binding
#         tag_00510.append(key)
#     elif genes_domains[key][0] == ['PF13632']: #Gtf2
#         tag_00520.append(key)
#     elif is_sublist(sorted(['PF17966', 'PF17965']), sorted(genes_domains[key][0])) or key == 'AKUH3B104X_PLPX00300': #LPXTG + mucin-binding
#         if ('PLPX' in key or int(key.replace('RS', '').split('_')[1]) < 1200 or 'MUB' in key) and key not in ['K2W83_RS00655', 'AAPFHON13_01040']:
#             tag_01020.append(key)
#             adh_dict[strain][0] += f'{key}, '
#         else:
#             tag_adh.append(key)
#             adh_dict[strain][1] += f'{key}, '
#         to_retrieve[strain].append(key)
#     elif genes_domains[key][0] == ['PF19087']: # and key.split('_')[0] not in ['K2W83', 'LDX55'] and int(key.split('_')[1]) < 1600:
#         tag_01250.append(key)
#     elif genes_domains[key][0] == ['PF19087']: #and key.split('_')[0] in ['K2W83', 'LDX55'] and int(key.split('_')[1].replace('RS', '')) < 800:
#         tag_01250.append(key)
#     elif sorted(genes_domains[key][0]) == sorted(['PF17966', 'PF00746']): #MubB2 + LPXTG
#         tag_14310.append(key)
    
# for strain in strains:
#     adh_dict[strain][0] = adh_dict[strain][0][:-2]
#     adh_dict[strain][1] = adh_dict[strain][1][:-2]
    
#     if adh_dict[strain][0] == '':
#         adh_dict[strain][0] = '-'
#     if adh_dict[strain][1] == '':
#         adh_dict[strain][1] = '-'

#     if adh_dict[strain] == ['-', '-']:
#         del adh_dict[strain]
    
#     if to_retrieve[strain] == []:
#         del to_retrieve[strain]
        

# with open(outfile, 'w') as handle:
#     handle.write('Strain\tAdhesin_locus1\tAdhesin_locus2\n')
#     for key in adh_dict.keys():
#         handle.write(f'{key}\t{adh_dict[key][0]}\t{adh_dict[key][1]}\n')
        
# with open(faa_outfile, 'w') as faa_out, open(faa_repset, 'w') as faa_rep:
#     for genome in to_retrieve.keys():
#         faa = f'{faa_dir}/{genome}_protein.faa'
#         check = False
#         if genome == 'HNS-8':
#             faa = faa.replace('/faa', '/new_genomes/faa')
#         if genome in repr_strains:
#             check = True
#         for record in SeqIO.parse(faa, 'fasta'):
#             if loctag_dict[record.id] in to_retrieve[genome]:
#                 record.id = loctag_dict[record.id]
#                 record.description = ''
#                 SeqIO.write(record, faa_out, 'fasta')
#                 if check:
#                     SeqIO.write(record, faa_rep, 'fasta')
            
        