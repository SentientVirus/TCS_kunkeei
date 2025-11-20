#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 17:06:38 2025

Script to get the locus tags of the proteins with YPKN or with LXET motifs in
the genome, plus the proteins with an LPXTG motif. It outputs: locus tag,
strain name, motif type (LPXTG/LXET or YPKN).

@author: Marina Mota-Merlo
"""

import os, re
from Bio import SeqIO
import multiprocessing
from functools import partial

def loctag2strain(locus_tag):
    prefix = locus_tag.split('_')[0].replace('AKU', '').replace('AAP', '')
    strain = prefix.replace('FHON', 'Fhon').replace('VQ058', 'GYUN-333')
    if strain.startswith('H'):
        strain = strain[:4] + '-' + strain[4:]
    strain = strain.replace('MUB42', 'HNS-8').replace('APS55', 'MP2')
    strain = strain.replace('K2W83', 'DSMZ').replace('LDX55', 'IBH001')
    return strain

class motif:
    '''
    Class to store the motif number, sequence and position.
    '''
    
    def __init__(self, mtype, seq, start, end): #Define motif class
        '''Requires:
            - mtype: motif type (str, YPKN/LXET/LPXTG)
            - seq: motif sequence (str, amino acids)
            - start: start position of the motif in the gene (int, in amino acids)
            - end: end position of the motif in the complete gene (int, in amino acids)'''
        self.mtype = mtype #Motif type
        self.seq = seq #Motif sequence
        self.start = start #Start position of the motif in a protein sequence
        self.end = end #End position of the motif in the sequence
        
    def __str__(self): #What print(class) returns
        return f'<Motif {self.mtype} from position {self.start} to {self.end}>'

class protein:
    '''
    Class to store protein objects with a series of motifs
    '''
    
    def __init__(self, locus_tag, length, motifs = []):
        
        '''
        Requires:
            - locus_tag: locus tag of the gene (str)
            - strain: name of the strain (str)
        Optional:
            - motifs: list of motif objects associated to the gene
            '''
        self.locus_tag = locus_tag
        self.length = length
        self.strain = loctag2strain(locus_tag)
        self.motifs = motifs
        
    def __str__(self):
        return f'<Locus {self.locus_tag} from strain {self.strain} with {len(self.motifs)} motifs of interest>'
    
    
def parse_files(file):
    prot_list = []
    loctag_list = []
    with open(file) as gbff:
        for record in SeqIO.parse(gbff, 'genbank'):
            features = [feature for feature in list(record.features) if 'translation' in feature.qualifiers.keys()]
            for feature in features:
                translation = feature.qualifiers['translation'][0]
                loctag = feature.qualifiers['locus_tag'][0]
                if loctag not in loctag_list:
                    new_prot = protein(loctag, len(translation), motifs = [])
                    loctag_list.append(loctag)
                for motifn in motif_list:
                    motif_info = re.finditer(motifn, translation)
                    info_list = list(motif_info)
                    if info_list != []:
                        for motif_match in list(info_list):
                            motseq = motif_match.group(0)
                            print(motseq)
                            start = motif_match.start() + 1 #Add one to the start (index starts at 0, but positions start at 1)
                            end = motif_match.end() + 1 #Add one to the end
                            mottype = motifn.replace('[A-Z]', 'x')
                            new_motif = motif(mottype, motseq, start, end)
                            new_prot.motifs.append(new_motif)
                if new_prot.motifs != []:
                    prot_list.append(new_prot)
    return prot_list

workdir = os.path.expanduser('~') + '/mucoid_project/adhesins'
indir = os.path.expanduser('~') + '/Akunkeei_files/gbff'
outfile = f'{workdir}/results/motifs/pilin_motifs.tsv'
threads = 4

if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))

infiles = [f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('.gbff') and 'M-0' not in file]
extra_files = [f'{indir.replace("gbff", "new_genomes/gbff")}/{file}' for file in os.listdir(indir.replace('gbff', 'new_genomes/gbff')) if file.endswith('.gbff')]
infiles += extra_files
# prot_dict = {}

motif_list = [r'YPKN', r'LP[A-Z]TG', r'Y[A-Z]L[A-Z]ET[A-Z]AP[A-Z]G']

if __name__ == '__main__': 
    pool = multiprocessing.Pool() 
    pool = multiprocessing.Pool(processes=threads)
    prot_list = pool.map(partial(parse_files), sorted(infiles))
    pool.map(parse_files, infiles)
    print("Input: {}".format(infiles))
    print("Output: {}".format(prot_list))
# prot_dict = parse_files(infiles)

# loctag_list = []
# prot_dict = {}

# for file in infiles:
#     with open(file) as gbff:
#         for record in SeqIO.parse(gbff, 'genbank'):
#             features = [feature for feature in list(record.features) if 'translation' in feature.qualifiers.keys()]
#             for feature in features:
#                 translation = feature.qualifiers['translation'][0]
#                 loctag = feature.qualifiers['locus_tag'][0]
#                 if loctag not in loctag_list:
#                     new_prot = protein(loctag, len(translation), motifs = [])
#                     loctag_list.append(loctag)
#                 for motifn in motif_list:
#                     motif_info = re.finditer(motifn, translation)
#                     info_list = list(motif_info)
#                     if info_list != []:
#                         for motif_match in list(info_list):
#                             motseq = motif_match.group(0)
#                             print(motseq)
#                             start = motif_match.start() + 1 #Add one to the start (index starts at 0, but positions start at 1)
#                             end = motif_match.end() + 1 #Add one to the end
#                             mottype = motifn.replace('[A-Z]', 'x')
#                             new_motif = motif(mottype, motseq, start, end)
#                             new_prot.motifs.append(new_motif)
#                 if new_prot.strain not in prot_dict.keys() and new_prot.motifs != []:
#                     prot_dict[new_prot.strain] = [new_prot]
#                 elif new_prot.motifs != []:
#                     prot_dict[new_prot.strain].append(new_prot)
    
with open(outfile, 'w') as motifs:
    motifs.write('Strain\tLocus_tag\tMotif_type\tMotif_sequence\tMotif_start\tMotif_end\tProtein_length\n')
    [motifs.write(f'{protein.strain}\t{protein.locus_tag}\t{motf.mtype}\t{motf.seq}\t{motf.start}\t{motf.end}\t{protein.length}\n') for sublist in prot_list for protein in sublist for motf in protein.motifs]
    # for sublist in prot_list:
    #     for 
    # for strain in prot_dict.keys():
    #     for protein in strain:
    #         for motf in protein:
    #             motifs.write(f'{strain}\t{protein.locus_tag}\t{motf.mtype}\t{motf.seq}\t{motf.start}\t{motf.end}\t{protein.length}\n')
