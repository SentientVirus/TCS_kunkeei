#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 17:11:34 2025

Script to retrieve the intergenetic regions between the MucBP adhesins, the
TF and other genes to search for a promoter.

@author: Marina Mota-Merlo
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import as_fasta

pos_tetR_tags = ['K2W83_RS00580', 'AKUFHON2_01070', #'AAPFHON13_00970', 
              'AKUG0101_01080', 'AKUG0102_01070', 'AKUG0103_01070', 
              'AKUG0401_01070', 'AKUG0402_01070', 'AKUG0403_PLPX00270',
              'AKUG0404_01070', 'AKUG0405_01070', 'AKUG0406_PLPX00280',
              'AKUG0407_01070', 'AKUG0408_01070', 'AKUG0410_01110',
              'AKUG0412_01110', 'AKUG0414_01070', 'AKUG0415_01070',
              'AKUG0417_01100', 'AKUG0420_PLPX00310', 'AKUG0601_01070',
              'AKUG0602_01070', 'AKUG0702_01070', 'AKUG0801_01070',
              'AKUG0802_01070', 'AKUG0803_01070', 'AKUG0804_01070',
              'AKUH1B104J_01070', 'AKUH1B105A_00980', 'AKUH3B101A_01060',
              'AKUH3B101J_01040', 'AKUH3B102A_01060', 'AKUH3B103J_01060',
              'AKUH3B103M_PLPX00270', 'AKUH3B104J_01040', 'AKUH3B104X_PLPX00270',
              'AKUH3B107A_01060', 'AKUH3B109M_01060', 'AKUH3B110M_01060',
              'AKUH3B111A_PLPX00260', 'AKUH3B111M_01060', 'AKUH3B202X_01050',
              'AKUH3B203J_01080', 'AKUH3B204J_01060', 'AKUH3B205J_01060',
              'AKUH3B207X_01060', 'AKUH3B208X_01070', 'AKUH4B202J_00980',
              'AKUH4B204J_01080', 'AKUH4B205J_01060', 'AKUH4B211M_01100',
              'AKUH4B412M_01140', 'AKUH4B501J_01140', 'AKUH4B502X_01070',
              'AKUH4B507J_01070', 'AKUH4B507X_01060', 'AKUH4B508X_01060',
              'MUB42_02675']

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

adh_tags = ['K2W83_RS00570', 'AKUFHON2_01050', #'AAPFHON13_00970', 
              'AKUG0101_01060', 'AKUG0102_01050', 'AKUG0103_01050', 
              'AKUG0401_01050', 'AKUG0402_01050', 'AKUG0403_PLPX00290',
              'AKUG0404_01050', 'AKUG0405_01050', 'AKUG0406_PLPX00300',
              'AKUG0407_01050', 'AKUG0408_01050', 'AKUG0410_01090',
              'AKUG0412_01090', 'AKUG0414_01050', 'AKUG0415_01050',
              'AKUG0417_01080', 'AKUG0420_PLPX00330', 'AKUG0601_01050',
              'AKUG0602_01050', 'AKUG0702_01050', 'AKUG0801_01050',
              'AKUG0802_01050', 'AKUG0803_01050', 'AKUG0804_01050',
              'AKUH1B104J_01050', 'AKUH1B105A_00960', 'AKUH3B101A_01040',
              'AKUH3B101J_01020', 'AKUH3B102A_01040', 'AKUH3B103J_01040',
              'AKUH3B103M_PLPX00290', 'AKUH3B104J_01020', 'AKUH3B104X_PLPX00290',
              'AKUH3B107A_01040', 'AKUH3B109M_01040', 'AKUH3B110M_01040',
              'AKUH3B111A_PLPX00280', 'AKUH3B111M_01040', 'AKUH3B202X_01030',
              'AKUH3B203J_01060', 'AKUH3B204J_01040', 'AKUH3B205J_01040',
              'AKUH3B207X_01040', 'AKUH3B208X_01050', 'AKUH4B202J_00960',
              'AKUH4B204J_01060', 'AKUH4B205J_01040', 'AKUH4B211M_01080',
              'AKUH4B412M_01120', 'AKUH4B501J_01120', 'AKUH4B502X_01050',
              'AKUH4B507J_01050', 'AKUH4B507X_01040', 'AKUH4B508X_01040',
              'MUB42_02660'] #1 loctag less for HNS-8

pre_adh_tags = ['K2W83_RS00565', 'AKUFHON2_01040', #'AAPFHON13_00970', 
              'AKUG0101_01050', 'AKUG0102_01040', 'AKUG0103_01040', 
              'AKUG0401_01040', 'AKUG0402_01040', 'AKUG0403_PLPX00300',
              'AKUG0404_01040', 'AKUG0405_01040', 'AKUG0406_PLPX00310',
              'AKUG0407_01040', 'AKUG0408_01040', 'AKUG0410_01080',
              'AKUG0412_01080', 'AKUG0414_01040', 'AKUG0415_01040',
              'AKUG0417_01070', 'AKUG0420_PLPX00340', 'AKUG0601_01040',
              'AKUG0602_01040', 'AKUG0702_01040', 'AKUG0801_01040',
              'AKUG0802_01040', 'AKUG0803_01040', 'AKUG0804_01040',
              'AKUH1B104J_01040', 'AKUH1B105A_00950', 'AKUH3B101A_01030',
              'AKUH3B101J_01010', 'AKUH3B102A_01030', 'AKUH3B103J_01030',
              'AKUH3B103M_PLPX00300', 'AKUH3B104J_01010', 'AKUH3B104X_PLPX00300',
              'AKUH3B107A_01030', 'AKUH3B109M_01030', 'AKUH3B110M_01030',
              'AKUH3B111A_PLPX00290', 'AKUH3B111M_01030', 'AKUH3B202X_01020',
              'AKUH3B203J_01050', 'AKUH3B204J_01030', 'AKUH3B205J_01030',
              'AKUH3B207X_01030', 'AKUH3B208X_01040', 'AKUH4B202J_00950',
              'AKUH4B204J_01050', 'AKUH4B205J_01030', 'AKUH4B211M_01070',
              'AKUH4B412M_01110', 'AKUH4B501J_01110', 'AKUH4B502X_01040',
              'AKUH4B507J_01040', 'AKUH4B507X_01030', 'AKUH4B508X_01030',
              'MUB42_02655'] #1 loctag less for HNS-8

PLPX_acc = ['OX335171.1', 'OX335123.1', 'OX335241.1', 'OX335188.1', 
            'OX335149.1', 'OX335191']

indir = os.path.expanduser('~') + '/Akunkeei_files/gbff'
indir2 = indir.replace('/gbff', '/new_genomes/gbff')
outfile = os.path.expanduser('~') + '/adhesins/sequences/intergenic/MucBP_intergenic.fna'
out_dir = os.path.dirname(outfile)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def get_unique_strains(tag_list):
    strains = []
    for tag in tag_list:
        strain = tag.split('_')[0].replace('AKU', '')
        if strain.startswith('H'):
            strain = strain[:4] + '-' + strain[4:]
        strain = strain.replace('K2W83', 'DSMZ12361').replace('MUB42', 'HNS-8').replace('FHON', 'Fhon')
        if strain not in strains:
            strains.append(strain)
    strains = sorted(strains)
    return strains

get_strains = get_unique_strains(adh_tags)

infiles = [f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('.gbff') and file.split('_')[0] in get_strains]
infiles += [f'{indir2}/{file}' for file in os.listdir(indir2) if file.endswith('.gbff') and file.split('_')[0] in get_strains]
infiles = sorted(infiles)
fna_infiles = [file.replace('gbff', 'fna') for file in infiles]

pos_dict = {}
for file in infiles:
    with open(file) as handle:
        for record in SeqIO.parse(handle, 'genbank'):
            for cds in record.features:
                if 'locus_tag' in list(cds.qualifiers.keys()):
                    loctag = cds.qualifiers['locus_tag'][0]
                    if loctag in pre_adh_tags and 'PLPX' not in loctag:
                        start = int(cds.location.end)
                    elif loctag in pre_adh_tags and 'PLPX' in loctag:
                        start = int(cds.location.start)
                    elif loctag in adh_tags and 'PLPX' not in loctag:
                        end = int(cds.location.start)
                    elif loctag in adh_tags and 'PLPX' in loctag:
                        end = int(cds.location.end)
    strain = os.path.basename(file).split('_')[0]
    pos_dict[strain] = (start, end)
    print(f'{strain}: {start}-{end}, {abs(end-start)} nucleotides.')
    
    
with open(outfile, 'w') as intergenic:
    for fna in fna_infiles:
        strain = os.path.basename(fna).split('_')[0]
        if pos_dict[strain][0] < pos_dict[strain][1]:
            plasmid = False
        else: plasmid = True
        with open(fna) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if not plasmid:
                    fna_seq = record.seq[pos_dict[strain][0]-30:pos_dict[strain][1]+30]
                    new_record = SeqRecord(fna_seq, id = f'{strain}_chromosome', description = '')
                    intergenic.write(as_fasta(new_record))
                    break
                elif 'plasmid: 2' in record.description:
                    fna_seq = record.seq[pos_dict[strain][1]-30:pos_dict[strain][0]+30].reverse_complement()
                    new_record = SeqRecord(fna_seq, id = f'{strain}_PLPX', description = '')
                    intergenic.write(as_fasta(new_record))
                    break
                
            
            