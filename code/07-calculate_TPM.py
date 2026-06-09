#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:12:54 2023

A script that calculates the number of TPM based on the number of counts

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Importing packages
# =============================================================================
import os, sys
import logging, traceback
import pandas as pd

# =============================================================================
# 0. Logging
# =============================================================================

logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
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

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# 1. Define input variables
# =============================================================================

outdir = os.path.dirname(snakemake.output[0]) #Output directory
outfile = snakemake.output.per_sample #Output file with per sample statistics
mean_outfile = snakemake.output.mean #Output file with per condition mean statistics
basedir = snakemake.params #Working directory

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

TPM = {} # Initialize dictionary for TPM calculation
    
# =============================================================================
# 2. Loop through the features in each of the files and calculate TPM
# =============================================================================
locus_tags = [] #List to store locus tags
for file in sorted(snakemake.input): #Loop through input files
    if len(locus_tags) == 0: #If the list with locus tags is empty
        check = True #Set the boolean to True
    else: check = False #Else, set it to False
    
    with open(file) as counts: #Open input file
        count_df = pd.read_csv(counts, sep = '\t', skiprows = 1) #Read the file as a dataframe
        T = 0 #Variable with the sum of normalized counts
        TPM_up = [] #List to store TPM scaled to the million
  
        for index, row in count_df.iterrows():
            rg = row.iloc[6] #Gene counts
            flg = row['Length']/10**3 #Gene length in kilobases
            TPM_calc = rg/flg #Gene counts normalized by gene length
            TPM_up.append(TPM_calc*1000000) #Normalized gene counts per million
            T += TPM_calc #Add normalized gene counts to the sum
            if check: #If the list of locus tags is empty
                locus_tags.append(row['Geneid']) #Add the locus tag of the gene

        basename = os.path.basename(file) #Get the file name
        TPM_list = [TPM_val/T for TPM_val in TPM_up] #Create a list with the final TPM values by dividing by the sum
        sample = f'I{basename[17:19]}_{basename[20:21]}_{basename[24:27].replace("_", "")}' #Get the name of the sample
        print(f'Calculated TPM for sample {sample}')
        TPM[sample] = TPM_list #Assign the TPMs to the sample

# =============================================================================
# 3. Save TPM calculations to files
# =============================================================================

df = pd.DataFrame.from_dict(TPM, orient='index').transpose() #Convert the TPM dictionary to a dataframe
df['locus_tag'] = locus_tags #Create a new column with locus tags
df = df.set_index('locus_tag') #Convert the new column into an index
df.to_csv(outfile, sep = '\t') #Save the dataframe to the output file as TSV
print('TPM saved to file') 

# =============================================================================
# 4. Calculate mean TPMs and save to file
# =============================================================================

isol_cond = [f'{i}{c}' for i in ['01', '02', '09', '10'] for c in ['F', 'S']] #Loop through isolates and media conditions
mean_df = pd.DataFrame() #Initialize an emoty dataframe
 
for ic in isol_cond: #Loop through isolates+condition
    underscore_ic = ic[:-1] + '_' +  ic[-1:] #Add an underscore between isolate and condition to match dataframe headers
    mean_df[ic] = df.filter(regex=underscore_ic).mean(axis=1) #Apply the mean function on columns with the same isolate+condition in the header
    mean_df.to_csv(mean_outfile, sep = '\t') #Write dataframe to file
print('Mean values saved to files')
