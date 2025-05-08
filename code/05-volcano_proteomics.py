#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 15:19:48 2025

@author: marina
"""

import os
import pandas as pd
from matplotlib import pyplot as plt
from math import log

workdir = os.path.expanduser('~') + '/proteomics'
indir = f'{workdir}/files/parsed'
outdir = f'{workdir}/plots'

if not os.path.exists(outdir):
    os.makedirs(outdir)

infiles = [f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('.tsv') and 'mucS_log' not in file]

for file in infiles:
    plt.figure(figsize= (15, 10))
    with open(file) as handle:
        df = pd.read_csv(handle, sep = '\t')
        pval = df.columns[-1]
        ratio = df.columns[-2]
        loctag = df.columns[1]
        label = df.columns[2]
        plot_title = ratio.replace('ratio_', '').replace('/', '_vs_').replace('_', ' ')
        filtered_df = df[df[ratio].notnull()]
        filtered_df = filtered_df[filtered_df[pval].notnull()]
        loctags = list(filtered_df[loctag])
        pvals = list(filtered_df[pval])
        ratios = list(filtered_df[ratio])
        log2fold = [log(fold, 2) for fold in ratios]
        log_pval = [-log(pv, 10) for pv in pvals]
        scatter_size = [8*(abs(fold)+2)/2 for fold in log2fold]
        colors = ['purple' if (abs(log2fold[i]) > 0.5 and log_pval[i] > 1) else 'grey' for i in range(len(log2fold))]
        labels =  list(df[label])
        final_labels = ['' if (labels[i] == '-' or abs(log2fold[i]) <= 0.5 or abs(log_pval[i]) <= 1) else labels[i] for i in range(len(log2fold))]
        print(file, log2fold, log_pval)
        plt.scatter(log2fold, log_pval, c = colors, s = scatter_size)
        plt.title(plot_title)
        plt.xlabel(r'$\mathregular{log_2 FC}$')
        plt.ylabel(r'$\mathregular{-log_10}$ p-value')
        for i in range(len(log2fold)):
            plt.annotate(final_labels[i], (log2fold[i], log_pval[i]))
            
        outfile = f'{outdir}/{plot_title.replace(" ", "_")}.png'
        plt.savefig(outfile)