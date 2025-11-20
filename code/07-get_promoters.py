#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 13:28:46 2025

@author: marina
"""

import os
import subprocess

workdir = os.path.expanduser('~') + '/mucoid_project/adhesins'
infile = f'{workdir}/sequences/intergenic/MucBP_intergenic.fna'
outdir = f'{workdir}/results/promoters'
threads = 12

subprocess.run(f'meme {infile} -oc {outdir} -dna -mod anr -nmotifs 20 -pal -w 6 -minsites 4', shell = True)
