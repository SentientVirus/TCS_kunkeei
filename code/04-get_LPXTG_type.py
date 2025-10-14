#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 18:29:09 2025

Script to check the type of LPXTG in each gene.

@author: marina
"""

import os
import re
from Bio import SeqIO

workdir = os.path.expanduser('~') + '/adhesins'
infile = f'{workdir}/sequences/Muc_adhesins.faa'

type_dict = {}
with open(infile) as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        LPXTG_motif = re.findall(r'LP(.{1,1})TG', str(record.seq))
        type_dict[record.id] = LPXTG_motif