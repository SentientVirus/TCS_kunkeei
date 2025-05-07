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

logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# 0. Define functions
# =============================================================================

def add_loctags(df, rep_dict):
    """Function to add a column with locus tags to a pre-existing dataframe
    Inputs: 
        - df: The input dataframe.
        - rep_dict: Dictionary with locus tags."""
        
    # df['Locus tags'] = df.loc[:, 'Majority protein IDs'] #Copy a column from the original dataframe
    df['Majority locus tags'] = df.loc[:, 'Majority protein IDs'] #Same as above
    
    for index, row in df.iterrows(): #Loop through the indexed rows of the dataframe
        # id_list = df.loc[index, 'Locus tags'].split(';') #This is needed if the Protein ID column can be kept
        # locus_list = [] #Create empty list to store locus tags
        # for protid in id_list: #Loop through proteins in the ID list
        
        for prot in replace_dict.keys(): #Loop through proteins in the dictionary
            rep = df.loc[index, 'Majority locus tags'] #Save the column to be modified
            if prot in rep: #If the protein ID is in the column
                df.loc[index, 'Majority locus tags'] = rep.replace(prot, replace_dict[prot]) #Replace it with the locus tag
                # if prot in protid: #If the protein ID is in the list of protein IDs
                #     locus = protid.replace(prot, replace_dict[prot]) #Replace the protein ID with the locus tag
        
            # locus_list.append(locus) #Add the current locus tag to the list
        # locus_tags = ';'.join(locus_list) #Join the list of locus tags into a single string
        # df.loc[index, 'Locus tags'] = locus_tags #Apply the changes to the column in the dataframe
        
    cols = list(df) #Get a list of the dataframe columns
    
    # cols.insert(0, cols.pop(cols.index('Locus tags')))
    cols.insert(0, cols.pop(cols.index('Majority locus tags'))) #Insert the locus tag column in the beginning
    
    df = df.loc[:, cols] #Apply changes to the dataframe
    
    return df #Return the modified dataframe

# =============================================================================
# 1. Define paths to inputs and outputs
# =============================================================================

workdir = os.path.expanduser('~') + '/proteomics' #Working directory
gbff = snakemake.input.gbk #os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff' #GenBank input with locus tag and protein ID information
# indir = f'{workdir}/files/parsed' #Directory with input files
signalP = snakemake.input.signalP #os.path.expanduser('~') + '/proteomics/results/SignalP/H3B1-04J_SignalP.txt' #Path to SignalP output
SP_tab = snakemake.output.signalP #signalP.replace('txt', 'tsv') #Path to SignalP tab file to be writted and loaded as a dataframe
outdir = os.path.dirname(snakemake.output.loci[0]) #f'{workdir}/files/loci' #Directory to save outputs
infiles = snakemake.input.infiles

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it
    
# =============================================================================
# 2. Retrieve information from the input files
# =============================================================================
    
replace_dict = {} #Create empty dictionary to store locus tags
with open(gbff) as handle: #Open GenBank file
    for record in GenBank.parse(handle): #Loop through records in the file
        for feature in record.features: #Loop through features in the records
            loctag = '' #Create a variable to store the locus tag
            prot_id = '' #Create a variable to store the protein ID
            for qual in feature.qualifiers: #Loop through qualifiers in the features
                if 'locus_tag' in qual.key: #If the qualifier is the locus tag
                    loctag = qual.value.strip('"') #Retrieve the locus tag removing any " at the end
                elif 'protein_id' in qual.key: #If the qualifier is the protein ID
                    prot_id = qual.value.strip('"') #Do the same as for the locus tag
            if loctag != '' and prot_id != '': #If there is a locus tag and a protein ID
                replace_dict[prot_id] = loctag #Store them as a key-value pair in a dictionary
                print(f'Saving protein {prot_id} with locus tag {loctag}...') #Print progress
                
                
with open(signalP) as SP_file, open(SP_tab, 'w') as tabfile: #Open the file with SignalP outputs
    for line in SP_file: #Loop through lines in the file
        line = line.split(' ') #Create a list by splitting the line by spaces
        while '' in line: #While the list contains empty strings
            line.remove('') #Remove the empty strings
        if 'name' in line: #If the list contains the string name
            line = line[1:] #Remove the first element of the list
            line[4] = 'end' #Update the value of the fifth element to the string end
        line_text = '\t'.join(line) + '\n' #Re-join the line as a tab-separated string, ending with a linebreak
        if '#' not in line_text: #If the character # is not in the string
            tabfile.write(line_text) #Write the string to a file
        
SP_df = pd.read_csv(SP_tab, sep = '\t') #Read the tab file as a dataframe
    
# =============================================================================
# 3. Add information to the input file and write it to an output
# =============================================================================

# infiles = [file for file in os.listdir(indir) if file.endswith('.tsv')] #Get a list of files to loop through
for infile in infiles: #Loop through the input files
    outfile = f'{outdir}/{os.path.basename(infile)}' #Set the path to the output file
    
    df = pd.read_csv(infile, sep = '\t') #Read the input file as a dataframe
    df = add_loctags(df, replace_dict) #Add column with the locus tags to the dataframe
    
    df['SP?'] = df.loc[:, 'Majority protein IDs'] #Copy a column from the original dataframe to create a signal peptide presence/absence column
    df['SP positions'] = df.loc[:, 'Majority protein IDs'] #Do the same to create a column to store signal peptide positions
    
    for index, row in df.iterrows(): #Loop through indexed rows in the dataframe from the MS input file
        for index2, row2 in SP_df.iterrows(): #Loop through rows in the SignalP dataframe
            if df.loc[index, 'Majority protein IDs'] == SP_df.loc[index2, 'name']: #If the locus tag in the MS dataframe matches the locus tag in SignalP
                print(f'Adding signal peptide prediction to protein {df.loc[index, "Majority protein IDs"]}') #Print progress
                df.loc[index, 'SP?'] = SP_df.loc[index2, '?'] #Add the SignalP information to the presence/absence column
                if SP_df.loc[index2, '?'] == 'Y': #If the signal peptide is present
                    df.loc[index, 'SP positions'] = f'1-{int(SP_df.loc[index2, "end"])-1}' #Add additional information to the position column
                else: #If it is not present
                    df.loc[index, 'SP positions'] = float('nan') #Add NaN to the position column
                
    cols = list(df) #Get a list of the dataframe columns
    cols.insert(4, cols.pop(cols.index('SP?'))) #Change the placement in the dataframe of the presence/absence column
    cols.insert(5, cols.pop(cols.index('SP positions'))) #Change the placement of the Signal peptide position column
    df = df.loc[:, cols] #Apply column changes to the dataframe
    
    df.to_csv(outfile, sep = '\t', index = False) #Save the dataframe to a tab-separated file
