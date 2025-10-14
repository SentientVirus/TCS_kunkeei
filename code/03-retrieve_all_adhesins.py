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

#Dictionary to color strains by phylogroup
phylogroup = {'A0901': 'C', 'A1001': 'F', 'A1002': 'C',
              'A1003': 'A', 'A1201': 'C', 'A1202': 'B',
              'A1401': 'B', 'A1404': 'E', 'A1802': 'C',
              'A1803': 'C', 'A1805': 'B', 'A2001': 'C',
              'A2002': 'F', 'A2003': 'F', 'A2101': 'B',
              'A2102': 'C', 'A2103': 'A', 'G0101': 'A', 
              'G0102': 'A', 'G0103': 'A', 'G0401': 'A',
              'G0402': 'A', 'G0403': 'B', 'G0404': 'A',
              'G0405': 'A', 'G0406': 'B', 'G0407': 'A',
              'G0408': 'A', 'G0410': 'A', 'G0412': 'A',
              'G0414': 'A', 'G0415': 'A', 'G0417': 'A',
              'G0420': 'B', 'G0601': 'A', 'G0602': 'A',
              'G0702': 'A', 'G0801': 'A', 'G0802': 'A',
              'G0803': 'A', 'G0804': 'A', 'Fhon2': 'A',
              'H1B1-04J': 'A', 'H1B1-05A': 'A', 'H1B3-02M': 'A',
              'H2B1-05J': 'C', 'H3B1-01A': 'A', 'H3B1-01J': 'A',
              'H3B1-01X': 'C', 'H3B1-02A': 'A', 'H3B1-02X': 'C',
              'H3B1-03J': 'A', 'H3B1-03X': 'C', 'H3B1-03M': 'A', 
              'H3B1-04J': 'A', 'H3B1-04X': 'A', 'H3B1-07A': 'A', 
              'H3B1-09M': 'A', 'H3B1-10M': 'A', 'H3B1-11A': 'A', 
              'H3B1-11M': 'A', 'H3B2-02M': 'C', 'H3B2-02X': 'A', 
              'H3B2-03J': 'A', 'H3B2-03M': 'C', 'H3B2-04J': 'A', 
              'H3B2-04M': 'C', 'H3B2-05J': 'A', 'H3B2-06M': 'C', 
              'H3B2-07X': 'A', 'H3B2-08X': 'A', 'H3B2-09X': 'B', 
              'H4B1-01A': 'C', 'H4B1-02A': 'C', 'H4B1-03J': 'C', 
              'H4B1-04A': 'C', 'H4B1-11J': 'C', 'H4B1-14J': 'C', 
              'H4B1-16J': 'C', 'H4B2-02J': 'A', 'H4B2-03M': 'C', 
              'H4B2-04J': 'A', 'H4B2-05J': 'A', 'H4B2-06J': 'C', 
              'H4B2-10M': 'C', 'H4B2-11M': 'A', 'H4B3-03J': 'C', 
              'H4B4-02J': 'A', 'H4B4-03J': 'C', 'H4B4-04J': 'C', 
              'H4B4-05J': 'C', 'H4B4-06M': 'C', 'H4B4-10M': 'C', 
              'H4B4-11M': 'C', 'H4B4-12M': 'A', 'H4B5-01J': 'A', 
              'H4B5-02X': 'A', 'H4B5-03X': 'A', 'H4B5-04J': 'B',
              'H4B5-05J': 'B', 'H4B5-07J': 'A', 'H4B5-07X': 'A', 
              'H4B5-08X': 'A', 'MP2': 'B', 'IBH001': 'C', 
              'DSMZ12361': 'A', 'HNS-8': 'X', 'GYUN-333': 'X', 'Fhon13': 'AAP'}

exclude = ['G0403', 'G0406', 'G0420', 'H3B1-02X', 'H3B2-03M', 'H4B1-02A', 
           'H4B4-05J', 'H4B4-10M']

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
            self.type = 'collagen-binding' #OK
            self.ref_locus = 'AKUH3B104J_00510'
        elif self.pfams == ['PF13632']:
            self.type = 'Gtf2' #OK
            self.ref_locus = 'AKUH3B104J_00520'
        elif self.pfams == ['PF19087', 'PF19087'] and 'LDX55' not in self.locus_tag and (self.locus_tag.split('_')[1].startswith('01') or self.locus_tag.split('_')[1].startswith('RS006') or self.locus_tag in ['MUB42_02735', 'APS55_RS07870']):
            self.type = 'SH3b' #OK
            self.ref_locus = 'AKUH3B104J_01250'
        elif is_sublist(sorted(['PF00746', 'PF19258']), sorted(self.pfams)):
            self.type = 'MubB2+LPXTG' #OK for now, but should be divided into subtypes
            self.ref_locus = 'AKUH3B104J_14310'
        elif is_sublist(sorted(['PF17966', 'PF17965']), self.pfams) or self.locus_tag == 'AKUH3B104X_PLPX00300':
            self.type = 'MucBP+LPXTG' #OK
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
outdir = f'{workdir}/results/adhesin_lists'
seqdir = f'{workdir}/sequences/adhesins'

pfam_domains = ['PF06458', 'PF19087', 'PF05737', 'PF19258', #MucBP, DUF5776, collagen-binding, signal peptide
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
                next_adh = adhesin(row[0], strain, pfams = [], domain_names = [])
                if next_adh.locus_tag != curr_adh.locus_tag:
                    curr_adh.get_type()
                    if curr_adh.locus_tag != '':
                        adh_dict[curr_adh.strain].append(curr_adh)
                    print(curr_adh)
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
        elif adhesin == 'MubB2+LPXTG':
            handle.write('Strain\tAdhesin_locus1\tAdhesin_locus2\tAdhesin_locus3\n')
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
            elif adhesin == 'MubB2+LPXTG':
                adhs = [adh for adh in adh_dict[strain] if adh.ref_locus == 'AKUH3B104J_14310']
                adhs.sort(key=lambda x: x.locus_tag)
                if len(adhs) == 2:
                    handle.write(f'{strain}\t{adhs[0].locus_tag}\t{adhs[1].locus_tag}\t-\n')
                elif len(adhs) > 0 and strain in exclude:
                    handle.write(f'{strain}\t{adhs[0].locus_tag}\t-\t-\n')
                elif len(adhs) > 0 and (phylogroup[strain] == 'F' or phylogroup[strain] == 'C'):
                    handle.write(f'{strain}\t-\t-\t{adhs[0].locus_tag}\n')
                elif len(adhs) > 0:
                    handle.write(f'{strain}\t-\t{adhs[0].locus_tag}\t-\n')
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
            
        