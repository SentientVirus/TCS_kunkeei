#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:41:24 2024

Script to add locus tags and SignalP predictions to the proteomics outputs

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================

import os, logging, traceback
import pandas as pd
from Bio import GenBank

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0]

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
# 0. Define functions
# =============================================================================

# def add_loctags(df, rep_dict):
#     """Function to add a column with locus tags to a pre-existing dataframe
#     Inputs: 
#         - df: The input dataframe.
#         - rep_dict: Dictionary with locus tags."""
        
#     df['Majority locus tags'] = df.loc[:, 'Majority protein IDs'] #Same as above
    
#     for index, row in df.iterrows(): #Loop through the indexed rows of the dataframe
        
#         for prot in replace_dict.keys(): #Loop through proteins in the dictionary
#             rep = df.loc[index, 'Majority locus tags'] #Save the column to be modified
#             if prot in rep: #If the protein ID is in the column
#                 df.loc[index, 'Majority locus tags'] = rep.replace(prot, replace_dict[prot]) #Replace it with the locus tag
        
#     cols = list(df) #Get a list of the dataframe columns
    
#     cols.insert(0, cols.pop(cols.index('Majority locus tags'))) #Insert the locus tag column in the beginning
    
#     df = df.loc[:, cols] #Apply changes to the dataframe
    
#     return df #Return the modified dataframe

# =============================================================================
# 1. Define paths to inputs and outputs
# =============================================================================

workdir = os.path.expanduser('~') + '/proteomics' #Working directory
signalP = snakemake.input.signalP #Path to SignalP output
SP_tab = snakemake.output.signalP #Path to SignalP tab file to be writted and loaded as a dataframe
outdir = os.path.dirname(snakemake.output.loci[0]) #Directory to save outputs
infiles = snakemake.input.infiles

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it
    
# =============================================================================
# 2. Retrieve information from the input files
# =============================================================================

with open(signalP) as SP_file, open(SP_tab, 'w') as tabfile: #Open the file with SignalP 
    SP_df = pd.read_csv(signalP, sep = '\t', skiprows = 1) #Read the file as a dataframe
    new_col = SP_df['# ID'].apply(lambda x: x.split(' ')[0]) #Create a new column only with protein IDs
    SP_df.pop('# ID') #Remove the previous ID column
    SP_df.insert(0, 'ID', new_col) #Add the new column to the dataframe at the first position
    SP_df.to_csv(SP_tab, sep = '\t', index = False) #Write the dataframe to a file without indexing it
    
# =============================================================================
# 3. Add information to the input file and write it to an output
# =============================================================================

# infiles = [file for file in os.listdir(indir) if file.endswith('.tsv')] #Get a list of files to loop through
for infile in infiles: #Loop through the input files
    outfile = f'{outdir}/{os.path.basename(infile)}' #Set the path to the output file
    
    df = pd.read_csv(infile, sep = '\t') #Read the input file as a dataframe
    
    df['SP prediction'] = df.loc[:, 'Majority protein IDs'] #Copy a column from the original dataframe to create a signal peptide presence/absence column
    df['SP positions'] = df.loc[:, 'Majority protein IDs'] #Do the same to create a column to store signal peptide positions
    
    for index, row in df.iterrows(): #Loop through indexed rows in the dataframe from the MS input file
        for index2, row2 in SP_df.iterrows(): #Loop through rows in the SignalP dataframe
            if df.loc[index, 'Majority protein IDs'] == SP_df.loc[index2, 'ID']: #If the locus tag in the MS dataframe matches the locus tag in SignalP
                print(f'Adding signal peptide prediction to protein {df.loc[index, "Majority protein IDs"]}') #Print progress
                df.loc[index, 'SP prediction'] = SP_df.loc[index2, 'Prediction'] #Add the SignalP information to the presence/absence column
                if SP_df.loc[index2, 'Prediction'] != 'OTHER': #If the signal peptide is present
                    df.loc[index, 'SP positions'] = f'1-{int(SP_df.loc[index2, "CS Position"].split(".")[0].split(" ")[2].split("-")[1])}' #Add additional information to the position column
                else: #If it is not present
                    df.loc[index, 'SP positions'] = float('nan') #Add NaN to the position column
                
    cols = list(df) #Get a list of the dataframe columns
    cols.insert(5, cols.pop(cols.index('SP prediction'))) #Change the placement in the dataframe of the presence/absence column
    cols.insert(6, cols.pop(cols.index('SP positions'))) #Change the placement of the Signal peptide position column
    df = df.loc[:, cols] #Apply column changes to the dataframe
    
    df.to_csv(outfile, sep = '\t', index = False) #Save the dataframe to a tab-separated file
