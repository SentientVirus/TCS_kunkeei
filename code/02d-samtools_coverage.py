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

#Create output directories if they don't exist        
[os.makedirs(newdir) for newdir in [outdir, plotdir] if not os.path.exists(newdir)]
    
name_dict = {'chromosome': 'OX335197.1', 'pKUN': 'OX335198.1'} #Dictionary with contig IDs
 
# =============================================================================
# 2. Calculate the sequencing depth and generate plots
# =============================================================================

for file in bamfiles: #Loop through BAM files
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
    outplot = f'{plotdir}/{isolate}_chromosome.png' #Path to output plot
    df_chr = df[df['contig'] == name_dict['chromosome']] #Retrieve data from the chromosome
    plt.plot('position', 'coverage', data = df_chr, color = '#7C55E6',
             linewidth = 1) #Plot the depth as a line
    plt.fill_between(df_chr['position'], df_chr['coverage'], #Fill the space under the line
                     where = df_chr['coverage'] >= 0, interpolate = False,
                     color = '#55BFE6')
    plt.margins(0, tight = True) #Remove plot margins
    plt.yticks([0, 100, 200, 300, 400, 500, 600, 700, 800]) #Set ticks of the y axis
    plt.title(f'Sequencing depth of isolate {isolate} (chromosome)') #Set plot title
    plt.xlabel('Position in the chromosome (bp)') #Set x axis label
    plt.ylabel('Read depth') #Set y axis label
    plt.savefig(outplot, dpi = 300) #Save to file with 300 dpi
    plt.savefig(outplot.replace('png', 'pdf'), dpi = 300) #Same, but in different formats
    plt.savefig(outplot.replace('png', 'svg'), dpi = 300)
    plt.show() #Show plot in the console
    
    outplot = f'{plotdir}/{isolate}_pKUN.png' #Same as above, but for the pKUN plasmid data
    df_pKUN = df[df['contig'] == name_dict['pKUN']]
    plt.plot('position', 'coverage', data = df_pKUN, color = '#E67C55', 
             linewidth = 1)
    plt.fill_between(df_pKUN['position'], df_pKUN['coverage'], 
                     where = df_pKUN['coverage'] >= 0, interpolate = False,
                     color = '#E6C355')
    plt.margins(0, tight = True)
    plt.yticks([0, 25, 50, 75, 100, 125, 150, 175, 200])
    plt.title(f'Sequencing depth of isolate {isolate} (pKUN)')
    plt.xlabel('Position in the plasmid (bp)')
    plt.ylabel('Read depth')
    plt.savefig(outplot, dpi = 300)
    plt.savefig(outplot.replace('png', 'pdf'), dpi = 300)
    plt.savefig(outplot.replace('png', 'svg'), dpi = 300)
    plt.show()