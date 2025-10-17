#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 16:47:34 2025

Script to retrieve the sequences of two genes neighboring the adhesins.

@author: Marina Mota-Merlo
"""

import os
import time
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#I added the tags for TetR manually, doesn't seem easy to retrieve from annotations

start_time = time.time()

inpath = os.path.expanduser('~') + '/Akunkeei_files/gbff'
workdir = os.path.expanduser('~') + '/adhesins'
outseqs = f'{workdir}/sequences/MucBP_neighbors'
outdir = outseqs.replace('sequences', 'alignments')
log = f'{workdir}/logs/03-get_neighboring_trees.log'
tree_dir = f'{workdir}/trees/MucBP_neighbors'
threads = 8

repr_strains = ['DSMZ12361', 'IBH001', 'GYUN-333', 'HNS-8', 'A0901', 
                'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805', 
                'Fhon2', 'G0102', 'G0403', 'H1B1-04J', 'H1B1-05A', 
                'H1B3-02M', 'H3B1-11M', 'H3B1-04J', 'H3B1-04X', 'H3B1-03M', #The last one is an extra strain, added because of the plasmid gene
                'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 'H4B2-06J', 'H4B4-02J', 
                'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 
                'H4B5-04J', 'H4B5-05J', 'MP2', 'Fhon13']

new_strains = ['GYUN-333', 'HNS-8']
extra_files = [inpath.replace('/gbff', '/new_genomes/gbff') + f'/{strain}_genomic.gbff' for strain in new_strains]
infiles = [f'{inpath}/{file}' for file in sorted(os.listdir(inpath)) if file.endswith('.gbff') and 'M-0' not in file]
infiles += extra_files
gene_dict = {}

genes = ['efpA', 'sasA', 'tetR', 'ipdC']
tag_dict = {'AAPFHON13_00810': 'sasA', 'AAPFHON13_01070': 'efpA', 'AAPFHON13_00980': 'ipdC',
            'K2W83_RS00530': 'sasA', 'K2W83_RS00670': 'efpA', 'K2W83_RS00580': 'ipdC',
            'APS55_RS02280': 'sasA', 'APS55_RS02155': 'efpA',
            'VQ058_RS00490': 'sasA', 'VQ058_RS00660': 'efpA',
            'MUB42_02625': 'sasA'} #efpA is labelled as a pseudogene in this strain

tetR_tags = ['K2W83_RS00575', 'AKUFHON2_01060', #'AAPFHON13_00970', 
              'AKUG0101_01070', 'AKUG0102_01060', 'AKUG0103_01060', 
              'AKUG0401_01060', 'AKUG0402_01060', 'AKUG0403_PLPX00280',
              'AKUG0404_01060', 'AKUG0405_01060', 'AKUG0406_PLPX00290',
              'AKUG0407_01060', 'AKUG0408_01060', 'AKUG0410_01100',
              'AKUG0412_01100', 'AKUG0414_01060', 'AKUG0415_01060',
              'AKUG0417_01090', 'AKUG0420_PLPX00320', 'AKUG0601_01060',
              'AKUG0602_01060', 'AKUG0702_01060', 'AKUG0801_01060',
              'AKUG0802_01060', 'AKUG0803_01060', 'AKUG0804_01060',
              'AKUH1B104J_01060', 'AKUH1B105A_00970', 'AKUH3B101A_01050',
              'AKUH3B101J_01030', 'AKUH3B102A_01050', 'AKUH3B103J_01050',
              'AKUH3B103M_PLPX00280', 'AKUH3B104J_01030', 'AKUH3B104X_PLPX00280',
              'AKUH3B107A_01050', 'AKUH3B109M_01050', 'AKUH3B110M_01050',
              'AKUH3B111A_PLPX00270', 'AKUH3B111M_01050', 'AKUH3B202X_01040',
              'AKUH3B203J_01070', 'AKUH3B204J_01050', 'AKUH3B205J_01050',
              'AKUH3B207X_01050', 'AKUH3B208X_01060', 'AKUH4B202J_00970',
              'AKUH4B204J_01070', 'AKUH4B205J_01050', 'AKUH4B211M_01090',
              'AKUH4B412M_01130', 'AKUH4B501J_01130', 'AKUH4B502X_01060',
              'AKUH4B507J_01060', 'AKUH4B507X_01050', 'AKUH4B508X_01050',
              'MUB42_02670']
        

[os.makedirs(out_dir) for out_dir in [outseqs, outdir, tree_dir] if not os.path.exists(out_dir)]

sasA_loctags = ['']
suffixes = ['', '_repset']
[open(f'{outseqs}/{gene}{suffix}.faa', 'w') for gene in genes for suffix in suffixes]

with open(log, 'w') as handle:
    handle.write('')

for file in infiles:
    strain = os.path.basename(file).split('_')[0]
    if strain in repr_strains:
        rep = True
    else: rep = False
    with open(file) as gbff:
        gbk = SeqIO.parse(gbff, 'genbank')
        for record in gbk:
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers.keys():
                    loctag = feature.qualifiers['locus_tag'][0]
                    if strain == 'H1B1-04J' and len(loctag.split('_')[1]) == 5 and 'R' not in loctag and 900 < int(loctag.split('_')[1]) < 1100:
                        print(loctag)
                    if (loctag in tetR_tags or 'gene' in feature.qualifiers.keys() or loctag in tag_dict.keys()) and 'translation' in feature.qualifiers.keys():
                        if loctag in tag_dict.keys():
                            gene_name = tag_dict[loctag]
                        elif loctag in tetR_tags:
                            gene_name = 'tetR'
                        else:
                            gene_name = feature.qualifiers['gene'][0]
                            
                        if gene_name == 'tetR' or ((gene_name == 'efpA' or gene_name == 'sasA' or gene_name == 'ipdC' or gene_name == 'kdc') and not (loctag.startswith('AKU') and int(loctag.split('_')[1]) > 2000)):
                            print(f'Strain: {strain}, locus: {loctag}, gene: {gene_name}')
                            gene_dict[loctag] = gene_name.replace('kdc', 'ipdC')
                            seq = feature.qualifiers['translation'][0]
                            new_record = SeqRecord(Seq(seq), id = loctag, 
                                                    name = gene_name,
                                                    description = '')
                            
                            with open(f'{outseqs}/{gene_name}{suffixes[0]}.faa', 'a') as all_faa:
                                SeqIO.write(new_record, all_faa, 'fasta')
                            if rep:
                                with open(f'{outseqs}/{gene_name}{suffixes[1]}.faa', 'a') as all_faa:
                                    SeqIO.write(new_record, all_faa, 'fasta')
                                
                            
for file in [f'{outseqs}/{gene}{suffix}.faa' for gene in genes for suffix in suffixes]:
    outfile = file.replace('.faa', '.mafft.faa').replace('sequences', 'alignments')
    subprocess.run(f'mafft-linsi --thread {threads} {file} > {outfile} 2>> {log};',
                    shell = True)
    subprocess.run(f'iqtree -nt AUTO -ntmax {threads} -redo -s {outfile} -st AA -msub nuclear -bb 1000 -bnni >> {log}', 
                    shell = True)
    subprocess.run(f'mv {outfile}.* {tree_dir}', shell = True)
    
end_time = time.time() - start_time
print(f'This script took {end_time/60:2f} minutes.')