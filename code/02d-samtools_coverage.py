#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 14:43:33 2026

Script to calculate per-position coverage using Samtools, and saving the results
to a file.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
import subprocess
import pandas as pd
from matplotlib import pyplot as plt

# =============================================================================
# 1. Set inputs and outputs
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/ugc00027'
indir = f'{workdir}/results/bam'
bamfiles = sorted([f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('.bam')]) #BAM files
outdir = f'{workdir}/results/coverage'
plotdir = f'{outdir}/plots'
bamlist = f'{indir}/bam_list.txt' #File including the path to the BAM files
outplot = f'{plotdir}/seqdepth.png'

#Create output directories if they don't exist        
[os.makedirs(newdir) for newdir in [outdir, plotdir] if not os.path.exists(newdir)]
    
name_dict = {'chromosome': 'OX335197.1', 'pKUN': 'OX335198.1'} #Dictionary with contig IDs

font_name = 'Arial'
titles = 14
plain = 12
 
# =============================================================================
# 2. Calculate the sequencing depth and generate plots
# =============================================================================

#Create the plot layout
fig, axs = plt.subplots(4, 2, sharex = 'col', sharey = 'col', 
                        figsize = (11.69, 8.27), dpi = 300) #A4, 300 DPI

fig.tight_layout(h_pad = 3, w_pad = 1) #Adjust spacing between subplots
fig.subplots_adjust(left = 0.1, top = 0.92) #Adjust figure margins

#Add a label for each isolate
fig.text(0.5, 0.95, 'Isolate 01', ha = 'center', va = 'center', 
         rotation = 'horizontal', weight = 'semibold', fontname = font_name,
         fontsize = titles)
fig.text(0.5, 0.715, 'Isolate 02', ha = 'center', va = 'center', 
         rotation = 'horizontal', weight = 'semibold', fontname = font_name,
         fontsize = titles)
fig.text(0.5, 0.478, 'Isolate 09', ha = 'center', va = 'center', 
         rotation = 'horizontal', weight = 'semibold', fontname = font_name,
         fontsize = titles)
fig.text(0.5, 0.245, 'Isolate 10', ha = 'center', va = 'center',
         rotation = 'horizontal', weight = 'semibold', fontname = font_name,
         fontsize = titles)

#Add a label for all the y axes (depth)
fig.text(0.05, 0.5, 'Sequencing depth', ha = 'center', va = 'center', 
         rotation = 'vertical', fontname = font_name, fontsize = plain)

for i in range(0, len(bamfiles)): #Loop through the list of BAM files
    file = bamfiles[i] #Retrieve the path to the file
    with open(bamlist, 'w') as handle: #Open the text file (overwrite mode)
        handle.write(f'{file}\n') #Write the BAM file
        
    isolate = os.path.basename(file).split('.')[0] #Retrieve isolate name
    print(f'Processing isolate {isolate}...')
    
    outfile = f'{outdir}/{isolate}.tab' #Output file
    
    #Run Samtools
    command = f'samtools depth -a -f {bamlist} -o {outfile} -q 0 -Q 10 -J -s' #Command to run
    subprocess.run(command, shell = True) #Run the command
    
    #Read results into a dataframe
    df = pd.read_csv(outfile, sep = '\t', header = None)
    df.columns = ['contig', 'position', 'coverage'] #Add column names
    
    #Plot results for the chromosome
    df_chr = df[df['contig'] == name_dict['chromosome']] #Retrieve data from the chromosome
    axs[i][0].plot('position', 'coverage', data = df_chr, color = '#7C55E6',
             linewidth = 1) #Plot the depth as a line
    axs[i][0].fill_between(df_chr['position'], df_chr['coverage'],
                     where = df_chr['coverage'] >= 0, interpolate = False,
                     color = '#55BFE6') #Fill the space under the line

    #Plot results for the plasmid (same as above)
    df_pKUN = df[df['contig'] == name_dict['pKUN']]
    axs[i][1].plot('position', 'coverage', data = df_pKUN, color = '#E67C55', 
             linewidth = 1)
    axs[i][1].fill_between(df_pKUN['position'], df_pKUN['coverage'], 
                     where = df_pKUN['coverage'] >= 0, interpolate = False,
                     color = '#E6C355')
    
    #Remove inner margins from both plots
    [axs[i][j].margins(0, tight = True) for j in range(0, len(name_dict.keys()))]

    if i == len(bamfiles)-1: #If it is the last row
        axs[i][0].set_xlabel('Position in the chromosome (bp)', #Add ax label
                             fontname = font_name, fontsize = plain)
        axs[i][1].set_xlabel('Position in the plasmid (bp)', #Same for pKUN
                             fontname = font_name, fontsize = plain)

plt.savefig(outplot, dpi = 300) #Save the plot to a figure (PNG)
plt.savefig(outplot.replace('png', 'pdf'), dpi = 300) #Save to PDF
plt.savefig(outplot.replace('png', 'svg'), dpi = 300) #Save to SVG