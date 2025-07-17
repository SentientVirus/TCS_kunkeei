#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 17:51:55 2025

Code to retrieve protein sequences from other A. kunkeei species from NCBI,
add the species name to the accession, and save the sequences to a single
fasta file.
Environment: ncbi_download.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import subprocess
import os
from Bio import SeqIO

# =============================================================================
# 1. Set path to output
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/blast_searches'
outfile = f'{workdir}/fasta/blastp_hits.faa'
finalfile = f'{workdir}/fasta/blastp_formatted.faa'

# =============================================================================
# 2. Define the sequences to be downloaded
# =============================================================================
acc_list = ['CAI2564285.1', 'WP_353317293.1', 'WP_120784742.1', 
            'WP_054658628.1', 'WP_220728881.1', 'WP_118903094.1', #Apilactobacillus and Bombilactobacillus
            'WP_089939307.1',  'WP_339959950.1', #Convivina and Nicoliella
            'WP_035452186.1', 'WP_034561863.1', 'WP_057737412.1',
            'HBH6538637.1', 'MTD40470.1', 'WP_349641026.1', 
            'WP_395391563.1', 'WP_010690715.1', #Fructobacillus
            'WP_407884581.1', 'WP_056974007.1', 'WP_123156317.1',
            'WP_280135435.1', 'WP_231533881.1', 'WP_125600976.1',
            'WP_413526436.1', 'WP_367304567.1', 'WP_048592711.1',
            'MCF7522510.1', 'WP_270333343.1', 'WP_031274475.1',
            'WP_224288405.1', 'WP_219480889.1', 'WP_013774648.1',
            'WP_002816652.1', 'WP_137602753.1', 'WP_058121442.1',
            'WP_205143364.1', 'WP_035441395.1', 'WP_168924558.1',
            'WP_014125860.1', 'WP_209527212.1', 'WP_063084080.1'
            ]

if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))
    
with open(outfile, 'w') as ncbi_out:
    ncbi_out.write('')

# =============================================================================
# 3. Create a string with the command to be run
# =============================================================================

list_command = ''
for acc in acc_list:
    list_command += f'{acc} '
list_command = list_command.strip()

full_command = f'ncbi-acc-download --molecule protein --out {outfile} {list_command}'

# =============================================================================
# 4. Run the command
# =============================================================================

subprocess.run(full_command, shell = True)

# =============================================================================
# 5. Read the output file and modify the fasta headers
# =============================================================================
    
records = []
with open(outfile) as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        species = record.description.split('[')[1]
        genus = species.split(' ')[0]
        epithet = species.split(' ')[1][:-1]
        if len(species.split(' ')) > 2 and epithet == 'sp':
            epithet += f'_{species.split(" ")[2]}_{species.split(" ")[3][:-1]}'
        elif len(species.split(' ')) > 2:
            epithet += f'i_{species.split(" ")[2][:-1]}'
        elif epithet == 'sp.':
            epithet = epithet[:-1] + '_G418'
        gtype = ''
        new_description = f'{record.id}_{gtype}{genus}_{epithet}'
        record.description = new_description
        record.id = new_description
        print(new_description)
        records.append(record)
        
with open(finalfile, 'w') as handle:
    SeqIO.write(records, handle, 'fasta')
        
