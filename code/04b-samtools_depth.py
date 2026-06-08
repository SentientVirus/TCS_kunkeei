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
import logging, traceback, sys
import subprocess
import pandas as pd
from matplotlib import pyplot as plt

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/snpseq00064'
log = f'{workdir}/logs/04b-plot_depth.log' #snakemake.log[0] #Path to log file

with open(log, 'w') as handle:
    handle.write('')

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    # Create a logger
    logger = logging.getLogger()

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Set inputs and outputs
# =============================================================================

bamdir = f'{workdir}/results/bam'
depth_dir = f'{workdir}/results/samtools_coverage'
bamfiles = sorted([f'{bamdir}/{file}' for file in os.listdir(bamdir) if file.endswith('.bam')]) #snakemake.input #Path to BAM files
indir = os.path.dirname(bamfiles[0]) #Input directory with the BAM files
#depthfiles = sorted([f'{depth_dir}/{file}' for file in os.listdir(depth_dir) if file.endswith('.cov')]) #snakemake.output.depth #Path to Samtools outputs
depthfiles = sorted([f'{depth_dir}/{file.replace(".bam", ".tsv")}' for file in os.listdir(bamdir) if file.endswith('.bam')])
outdir = os.path.dirname(depthfiles[0]) #Directory with the Samtools outputs
outplots = [f'{workdir}/plots/coverage/coverage_plot.png', #snakemake.output.plots #Output plot
            f'{workdir}/plots/coverage/coverage_plot.pdf',
            f'{workdir}/plots/coverage/coverage_plot.svg']
plotdir = os.path.dirname(outplots[0]) #Directory with the output plot
bamlist = f'{indir}/bam_list.txt' #File including the path to the BAM files

#Create output directories if they don't exist        
[os.makedirs(newdir) for newdir in [outdir, plotdir] if not os.path.exists(newdir)]
    
name_dict = {'chromosome': 'OX335197.1', 'pKUN': 'OX335198.1'} #Dictionary with contig IDs

#Text settings
font_name = 'Arial' #Font type
titles = 14 #Font size (titles)
plain = 12 #Font size (body)
no_replicates = 5 #Number of replicates

plot_conditions = {'Mucoid +S': 0, 'Mucoid -S': 1, 'Aggregating +S': 2,
                   'Aggregating -S': 3}

isolate_colors = {'01': 'magenta', '02': 'orange', '09': 'purple',
                  '10': 'cyan'}
 
# =============================================================================
# 2. Calculate the sequencing depth and generate plots
# =============================================================================

print('Formatting subplots...')
#Create the plot layout
fig, axs = plt.subplots(4, 2, sharex = 'col', #sharey = 'col', 
                        figsize = (11.69, 8.27), dpi = 300) #A4, 300 DPI

fig.tight_layout(h_pad = 2, w_pad = 1) #Adjust spacing between subplots
fig.subplots_adjust(left = 0.1, top = 0.92, bottom = 0.1) #Adjust figure margins

#Add a label for all the y axes (depth)
fig.text(0.05, 0.5, 'Sequencing depth', ha = 'center', va = 'center', 
         rotation = 'vertical', fontname = font_name, fontsize = plain)

for i in range(0, len(bamfiles)): #Loop through the list of BAM files
    file = bamfiles[i] #Retrieve the path to the file
    with open(bamlist, 'w') as handle: #Open the text file (overwrite mode)
        handle.write(f'{file}\n') #Write the BAM file
        
    sample = os.path.basename(file).split('.')[0] #Retrieve isolate name
    print(f'Processing sample {sample}... ({i+1}/{len(bamfiles)}')
    
    isolate = sample.split('-')[4] #Retrieve isolate number from file name
    condition = sample.split('-')[5] #Retrieve condition from file name
    
    #Assign a phenotype based on isolate
    if isolate == '01' or isolate == '02':
        phenotype = 'Mucoid'
    else: phenotype = 'Aggregating'
    
    #Assign a condition based on file name
    if condition == 'S':
        condition = '+S'
    elif condition == 'F':
        condition = '-S'
        
    sample_cond = f'{phenotype} {condition}' #Create a variable with sample and condition information
    
    outfile = depthfiles[i] #Output file
    
    #Run Samtools
    command = f'samtools depth -a -f {bamlist} -o {outfile} -q 0 -Q 10 -J -s' #Command to run
    subprocess.run(command, shell = True) #Run the command
    
    print(f'Samtools results saved to {outfile}!')
    
    #Read results into a dataframe
    df = pd.read_csv(outfile, sep = '\t', header = None)
    df.columns = ['contig', 'position', 'coverage'] #Add column names
    
    #Plot results for the chromosome
    print(f'Plotting the depth of the chromosome of isolate {isolate}...')
    df_chr = df[df['contig'] == name_dict['chromosome']] #Retrieve data from the chromosome
    axs[plot_conditions[sample_cond]][0].plot('position', 'coverage', data = df_chr, 
                                              color = isolate_colors[isolate], #'#7C55E6',
                                              alpha = 0.1, linewidth = 1) #Plot the depth as a line
    axs[plot_conditions[sample_cond]][0].set_title(f'{sample_cond}, chromosome')

    #Plot results for the plasmid (same as above)
    print(f'Plotting the depth of the pKUN of isolate {isolate}...')
    df_pKUN = df[df['contig'] == name_dict['pKUN']]
    axs[plot_conditions[sample_cond]][1].plot('position', 'coverage', data = df_pKUN, 
                                              color =  isolate_colors[isolate], #'#E67C55',
                                              alpha = 0.1, linewidth = 1)
    axs[plot_conditions[sample_cond]][1].set_title(f'{sample_cond}, pKUN')
    
    #Remove inner margins from both plots
    [axs[plot_conditions[sample_cond]][j].margins(0, tight = True) for j in range(0, len(name_dict.keys()))]

    if plot_conditions[sample_cond] == max(plot_conditions.values()): #If it is the last row
        axs[plot_conditions[sample_cond]][0].set_xlabel('Position in the chromosome (Mbp)', #Add ax label
                             fontname = font_name, fontsize = plain)
        axs[plot_conditions[sample_cond]][1].set_xlabel('Position in the plasmid (bp)', #Same for pKUN
                             fontname = font_name, fontsize = plain)

print('Removing intermediate file...')
subprocess.run(f'rm {bamlist}', shell = True) #Remove the file with the path to the BAM

print('Saving plot...')
#Save the plot to a figure (PNG, PDF and SVG)
[plt.savefig(outplot, dpi = 300) for outplot in outplots]
print('Done!')