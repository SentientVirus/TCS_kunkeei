#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 17:17:10 2026

Script to parse the Perseus outputs and make them easier to process in R.

This script is dependent on steps 01 and 02 of the main proteomics pipeline,
and step 03a is running Perseus separately.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
import pandas as pd
import logging, sys

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/proteomics'
log = f'{workdir}/logs/04a-parse_Perseus.py'

with open(log, 'w') as logfile: #Overwrite log file
    logfile.write('')
            
#Redirect stdout and stderr to log file
sys.stdout = open(log, 'a')
sys.stderr = open(log, 'a')

#Format the logging
logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================

infolder = f'{workdir}/Perseus/results' #Input directory
outfolder = f'{infolder}/parsed' #Output directory

#Create output directory if it doesn't exist
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

#List of input files
infiles = [f'{infolder}/{file}' for file in os.listdir(infolder) if file.endswith('.txt') and 'mucS_log' not in file]
#List of output files
outfiles = [f'{outfolder}/{file.replace(".txt", ".tsv")}' for file in os.listdir(infolder) if file.endswith('.txt') and 'mucS_log' not in file]

# =============================================================================
# 2. Loop through files and create outpus
# =============================================================================

for i in range(len(infiles)): #Loop through the length of the input files
    #Retrieve the name of the...
    file = infiles[i] #...input file
    outfile = outfiles[i] #...output file
    
    #Name of the comparison
    comparison = os.path.basename(file).replace('.txt', '')
    rev = comparison.split('_') #Split the name
    
    #Retrieve the comparison name, as it appears in the column names
    rev_comparison = f'{rev[3]}_{rev[4]}_{rev[0]}_{rev[1]}'
    
    df = pd.read_csv(file, sep = '\t') #Read the input file as a dataframe
    
    #Get the dataframe columns
    col_list = df.columns
    
    #Create a dictionary to use the desired column names
    new_cols = {col: col.split(':')[-1].replace(rev_comparison, '').replace('significant', 'comparison').strip() for col in col_list}
    df.rename(columns = new_cols, inplace = True) #Replace column names
    
    #Set the desired column order
    col_order = ['Majority locus tags', 'Majority protein IDs',
                 'Original annotations', 'Refined annotations',
                 'SP prediction', "Welch's T-test Significant",
                 "Welch's T-test comparison", 'Number Of Imputations',
                 "-Log Welch's T-test p-value", "Welch's T-test q-value",
                 "Welch's T-test Difference", 
                 "Welch's T-test Test statistic"] + list(col_list[:6])
    
    #Create a reordered dataframe
    df_ordered = df[col_order]
    
    #Save the dataframe to file
    df_ordered.to_csv(outfile, sep = '\t', index = False)