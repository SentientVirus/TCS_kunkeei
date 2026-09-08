#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 16:22:43 2024

Script to compare the results of Phobius to those of DeepTMHMM and keep
the overlapping transmembrane helices, then save them to a file together
with the protein ID and protein locus tag.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries
# =============================================================================

import os
import pandas as pd
from Bio import GenBank
import logging, sys

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
# 1. Create a class of objects to store TM prediction information
# =============================================================================

class TM:
    def __init__(self, prot_id, locus_tag, start, end):
        self.prot_id = prot_id #Protein ID from NCBI
        self.locus_tag = locus_tag #Locus tag of the gene
        self.start = int(start) #Start position of the gene
        self.end = int(end) #End position of the gene
        self.length = self.end - self.start #Gene length
    def __str__(self): #Function to print the object
        return f'TMH in protein {self.prot_id} ({self.start}-{self.end})' #Information to be printed
    def write_TM(self, out_file): #Function to write the object to a file
        with open(out_file, 'a') as handle: #Open the target file in append mode
            handle.write(f'{self.locus_tag}\t{self.prot_id}\t{self.start}\t{self.end}\n') #Add TMH information in a tab-separated line

# =============================================================================
# 2. Define input variables and paths
# =============================================================================

gbff = snakemake.input.gbk #os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff' #Reference GenBank file
phob = snakemake.input.Phobius #f'{workdir}/results/Phobius/H3B1-04J_phobius.txt' #Path to Phobius predictions
DTMH = snakemake.input.DeepTMHMM #f'{workdir}/results/DeepTMHMM/TMRs.gff3' #Path to DeepTMHMM predictions
pred_out = snakemake.output.common_pred #f'{workdir}/results/TMH_predictions/TMH.tab' #File to save predictions that are consistent between methods 
infiles = snakemake.input.sample_in #Input files with proteomics results and SP predictions
outfiles = snakemake.output.sample_pred #Output files

if not os.path.exists(os.path.dirname(pred_out)): #If the directory where the consistent predictions will be saved does not exist
    os.makedirs(os.path.dirname(pred_out)) #Create it
 
# =============================================================================
# 3. Open the input files and retrieve all the information
# =============================================================================

with open(gbff) as handle: #Open the input GenBank file
    logging.info('Loading GenBank file...')
    replace_dict = {} #Create an empty dictionary
    for record in GenBank.parse(handle): #Loop through records in the GenBank file
        for feature in record.features: #Loop through features in the GenBank file
            loctag = '' #Empty locus tag
            prot_id = '' #Empty protein ID
            for qual in feature.qualifiers: #Loop through feature qualifiers
                if 'locus_tag' in qual.key: #If the qualifier is the locus tag
                    loctag = qual.value.strip('"') #Retrieve the locus tag, removing " at the end
                elif 'protein_id' in qual.key: #If the qualifier is the protein ID
                    prot_id = qual.value.strip('"') #Retrieve the protein ID, removing " at the end
            if loctag != '' and prot_id != '': #If there's a locus tag and a protein ID
                replace_dict[prot_id] = loctag #Add the locus tag to the dictionary
                logging.info(f'Saving protein {prot_id} with locus tag {loctag}...') #Print progress
    
with open(phob) as phobius: #Open Phobius output
    logging.info('Loading Phobius results...') 
    phobius_TM = [] #List to store Phobius predictions
    for line in phobius: #Loop through Phobius predictions
        if line.startswith('ID'): #If the line includes an ID
           line = line.strip() #Remove special characters at the end
           prot_id = line.split(' ')[-1] #Divide the line and extract the protein ID
           loctag = replace_dict[prot_id] #Retrieve the locus tag
               
        elif 'TRANSMEM' in line: #If the prediction is a transmembrane helix
            line = line.strip() #Remove special characters at the end of the line
            pos = line.split(' ') #Split the line into a list by spaces
            while '' in pos: #If there are empty strings in the list
                pos.remove('') #Remove empty strings from the list
            helix = TM(prot_id, loctag, pos[2], pos[3]) #Create TM object
            logging.info(helix) #Print object
            phobius_TM.append(helix) #Add object to the Phobius list
    
with open(DTMH) as deepTM: #Open DeepTMHMM output
    logging.info('Loading DeepTMHMM results...')
    deep_TM = [] #List to store DeepTMHMM predictions
    for line in deepTM: #Loop through DeepTMHMM predictions
        if line.startswith('CAI') and 'TMhelix' in line: #If the prediction is a transmembrane helix
            line = line.strip() #Remove special characters at the end of the line
            line = line.split('\t') #Split line into a list by tabs
            while '' in line: #If there are empty strings in the list
                line.remove('') #Remove empty strings from the list
            helix = TM(line[0], replace_dict[line[0]], line[2], line[3]) #Create helix object
            logging.info(helix) #Print object
            deep_TM.append(helix) #Add object to the Phobius list
          
logging.info('Save results to dictionary...')
add_info = {} #Dictionary with the information to be incorporated into the output file
for TMH1 in phobius_TM: #Loop through helices in Phobius predictions
    for TMH2 in deep_TM: #Loop through helices in DeepTMHMM predictions
        if TMH1.prot_id == TMH2.prot_id and (TMH1.start < TMH2.end and TMH1.end > TMH2.start): #If the predictions overlap
            new_helix = TM(TMH1.prot_id, TMH1.locus_tag, max(TMH1.start, TMH2.start), min(TMH1.end, TMH2.end)) #Save the overlap as a new helix
            new_helix.write_TM(pred_out) #Write the new helix to the output prediction file
            if new_helix.prot_id not in add_info.keys(): #If the protein with helices hasn't been incorporated into the dictionary
                add_info[new_helix.prot_id] = f'{new_helix.start}-{new_helix.end}' #Incorporate it with the helix information
            else: #Otherwise
                add_info[new_helix.prot_id] += f';{new_helix.start}-{new_helix.end}' #Append the helix information to preexisting information
    
# =============================================================================
# 4. Add TMH predictions
# =============================================================================

logging.info('Add predictions...')
for i in range(len(infiles)): #Loop through input files
    outfile = outfiles[i] #Set the path and name of the output file
    
    logging.info(f'Input file: {infiles[i]}\nOutput file: {outfile}')
    
    with open(outfile, 'w') as out_pred: #Open the output file in write mode
        out_pred.write('locus_tag\tprotein_id\tstart\tend\n') #Write the file headers
    
    df = pd.read_csv(infiles[i], sep = '\t') #Read the input file
    
    df['TMH?'] = df.loc[:, 'Majority protein IDs'] #Create a new column in the dataframe by copying another column
    df['TMH positions'] = df.loc[:, 'Majority protein IDs'] #Create a new column in the dataframe by copying another column
    for index, row in df.iterrows(): #Loopthrough index and row in the dataframe
        if df.loc[index, 'Majority protein IDs'] in add_info.keys(): #If the majority protein ID in the row is in the keys of the helix dictionary
            df.loc[index, 'TMH?'] = 'Y' #Set the value of the presence/absence column to Y(es)
            df.loc[index, 'TMH positions'] = add_info[df.loc[index, 'Majority protein IDs']] #Set the value of the column with helix positions to 
        else: #Otherwise
            df.loc[index, 'TMH?'] = 'N' #Set the value of the presence/absence column to N(o)
            df.loc[index, 'TMH positions'] = float('nan') #Set the positions to NaN
            
    cols = list(df) #Get the dataframe columns
    
    cols.insert(7, cols.pop(cols.index('TMH?'))) #Change the position of the new columns
    cols.insert(8, cols.pop(cols.index('TMH positions')))
    df = df.loc[:, cols] #Apply changes to the dataframe
    
    df.to_csv(outfile, sep = '\t', index = False) #Save the dataframe to a tab-separated file´
        



