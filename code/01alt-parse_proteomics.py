#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 18:28:49 2026

Script that reads all the Excel files with proteomics data, converts them into
tables and generates a file with all the LFQ values of each sample
       
The output file is the input of other scripts to add locus tags, signal peptide
and TM helices predictions.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries
# =============================================================================

import numpy as np
import pandas as pd
import os, logging, sys
from Bio import SeqIO

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/proteomics'

log = f'{workdir}/logs/01alt-parse_proteomics.log'

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
# 1. Define paths to inputs and other variables
# =============================================================================

outdir = f'{workdir}/results/parsed_MS' #snakemake.params.outdir #Output directory

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

exp1 = f'{workdir}/files/MS-24-030_MaxQuant_results.xlsx' #snakemake.input.exp1 #File with the results of the first MS dataset

exp2 = f'{workdir}/files/MS-24-038_MaxQuant_results.xlsx' #snakemake.input.exp2 #File with the results of the second MS dataset

#Original annotations (annot1) and refined annotations (annot2)
annot1 = os.path.expanduser('~') + '/Akunkeei_files/gbff/H3B1-04J_genomic.gbff'
annot2 = os.path.expanduser('~') + '/mucoid_project/ugc00027/results/annotations/emapper2gbk/reference.gbk' #snakemake.input.gbk #Reference GenBank file

exp1_dict = {1: 'mucSD_stat01', 2: 'mucSD_stat02', 3: 'mucSD_stat03', #Sample annotations of the first MS dataset
             4: 'inhS_stat04', 5: 'inhS_stat05', 6: 'inhS_stat06'}

exp2_dict = {1: 'mucS_log01', 2: 'mucS_log02', 3: 'mucS_log03', #Sample annotations of the second dataset
             4: 'mucF_log04', 5: 'mucF_log05', 6: 'mucF_log06',  
             7: 'mucF_stat07', 8: 'mucF_stat08', 9: 'mucF_stat09',
             10: 'mucFD_stat10', 11: 'mucFD_stat11', 12: 'mucFD_stat12',
             13: 'inhS_log13', 14: 'inhS_log14', 15: 'inhS_log15',
             16: 'inhF_log16', 17: 'inhF_log17', 18: 'inhF_log18',
             19: 'inhF_stat19', 20: 'inhF_stat20', 21: 'inhF_stat21'}

# manual_annot = {'CAI2650094.1': 'GS-BrS', 'CAI2649694.1': 'GS', 
#                 'CAI2649852.1': 'BrS', 'CAI2649985.1': 'BrS',
#                 'CAI2554085.1': 'adhesin_1020', 'CAI2550733.1': 'gtf2_520',
#                 'CAI2550731.1': 'adhesin_510', 'CAI2659529.1': 'adhesin_13900',
#                 'CAI2663438.1': 'adhesin_14310', 'CAI2671120.1': 'parA',
#                 'CAI2671199.1': 'repC', 'CAI2671279.1': 'kukP',
#                 'CAI2671341.1': 'kukA', 'CAI2671425.1': 'kukC',
#                 'CAI2671510.1': 'kukT', 'CAI2671601.1': 'kukF',
#                 'CAI2671711.1': 'kukE', 'CAI2671757.1': 'kukG1',
#                 'CAI2671762.1': 'kukG2', 'CAI2671872.1': 'kukB',
#                 'CAI2671986.1': 'pKUN_HK', 'CAI2672095.1': 'HTH',
#                 'CAI2672155.1': 'tnpR', 'CAI2672555.1': 'repA',
#                 'CAI2555830.1': '01250', 'CAI2561695.1': 'giant1',
#                 'CAI2561700.1': 'giant2', 'CAI2561781.1': 'giant3',
#                 'CAI2561859.1': 'giant4', 'CAI2561935.1': 'giant5'}

cds_dict1 = {} #Create an empty dictionary to store the original annotations
cds_dict2 = {} #Same, but for the modified annotations
prot2loctag = {} #Dictionary to assign locus tags to the protein IDs

# =============================================================================
# 2. Read input data
# =============================================================================

logging.info(f'Loading {exp1} dataset...')
exp1_df = pd.read_excel(exp1, sheet_name = 'All') #Read results from MS-24-030 experiment
exp1_df = exp1_df.drop('i', axis = 1) #Remove column i
exp1_df = exp1_df[exp1_df['Majority protein IDs'].notna()] #Filter out empty rows and rows with column averages
max_index = max(exp1_dict.keys()) + 1 #Get the blank column index
exp1_df[f'LFQ intensity {max_index}'] = exp1_df[f'LFQ intensity {max_index}'].fillna(0) #Change the nan of the blank to 0s (otherwise, n - nan = nan)
exp1_df[f'Peptides {max_index}'] = exp1_df[f'Peptides {max_index}'].fillna(0)

logging.info(f'Loading {exp2} dataset...')
exp2_df = pd.read_excel(exp2, sheet_name = 'All')
[exp2_df.insert(25+i, f'LFQ intensity {i}', np.nan) for i in range(1, 3)] #Insert missing samples as NaNs
[exp2_df.insert(5+i, f'Peptides {i}', np.nan) for i in range(1, 3)] #Insert missing samples as NaNs
exp2_df = exp2_df[exp2_df['Majority protein IDs'].notna()]
max_index = max(exp2_dict.keys()) + 1
exp2_df[f'LFQ intensity {max_index}'] = exp2_df[f'LFQ intensity {max_index}'].fillna(0) #Change the nan of the blank to 0s (otherwise, n - nan = nan)
exp2_df[f'Peptides {max_index}'] = exp2_df[f'Peptides {max_index}'].fillna(0)

# =============================================================================
# 3. Create two dictionaries with annotations
# =============================================================================

logging.info(f'Loading original annotations from {annot1}...')
with open(annot1) as handle: #Open annotation file
    gbk = SeqIO.parse(handle, 'genbank') #Parse file as GenBank
    for record in gbk: #Loop through records in the file
        for cds in record.features: #Loop through features in each record
            if 'locus_tag' in cds.qualifiers.keys() and 'protein_id' in cds.qualifiers.keys():
                prot2loctag[cds.qualifiers['protein_id'][0]] = cds.qualifiers['locus_tag'][0]
            if cds.type == 'CDS' and 'gene' in cds.qualifiers: #If the feature is a CDS and there is a gene name in the annotations
                cds_dict1[cds.qualifiers['protein_id'][0]] = cds.qualifiers['gene'][0] #Add the gene name to the dictionary, using the protein ID as a key
            elif cds.type == 'CDS' and cds.qualifiers['protein_id'][0] not in cds_dict1.keys(): #If the gene name is not in the annotations
                cds_dict1[cds.qualifiers['protein_id'][0]] = '-' #Use a minus symbol as gene name

logging.info(f'Loading refined annotations from {annot2}...')
with open(annot2) as handle:
    gbk = SeqIO.parse(handle, 'genbank')
    for record in gbk:
        for cds in record.features:
            if cds.type == 'CDS' and 'gene' in cds.qualifiers: #If the feature is a CDS and there is a gene name in the annotations
                cds_dict2[cds.qualifiers['locus_tag'][0]] = cds.qualifiers['gene'][0] #Add the gene name to the dictionary, using the protein ID as a key
            elif cds.type == 'CDS' and cds.qualifiers['locus_tag'][0] not in cds_dict2.keys(): #If the gene name is not in the annotations
                cds_dict2[cds.qualifiers['locus_tag'][0]] = '-' #Use a minus symbol as gene name
                
print(cds_dict2)
                
logging.info('Adding manual annotations to the refined annotations...')
#Manually add the annotations of proteins of interest
cds_dict2['CAI2650094.1'] = 'GS-BrS'
cds_dict2['CAI2649694.1'] = 'GS'
cds_dict2['CAI2649852.1'] = 'BrS'
cds_dict2['CAI2649985.1'] = 'BrS'
cds_dict2['CAI2554085.1'] = 'adhesin_1020'
cds_dict2['CAI2550733.1'] = 'gtf2_520'
cds_dict2['CAI2550731.1'] = 'adhesin_510'
cds_dict2['CAI2659529.1'] = 'adhesin_13900'
cds_dict2['CAI2663438.1'] = 'adhesin_14310'
cds_dict2['CAI2671120.1'] = 'parA'
cds_dict2['CAI2671199.1'] = 'repC'
cds_dict2['CAI2671279.1'] = 'kukP'
cds_dict2['CAI2671341.1'] = 'kukA'
cds_dict2['CAI2671425.1'] = 'kukC'
cds_dict2['CAI2671510.1'] = 'kukT'
cds_dict2['CAI2671601.1'] = 'kukF'
cds_dict2['CAI2671711.1'] = 'kukE'
cds_dict2['CAI2671757.1'] = 'kukG1'
cds_dict2['CAI2671762.1'] = 'kukG2'
cds_dict2['CAI2671872.1'] = 'kukB'
cds_dict2['CAI2671986.1'] = 'pKUN_HK'
cds_dict2['CAI2672095.1'] = 'HTH'
cds_dict2['CAI2672155.1'] = 'tnpR'
cds_dict2['CAI2672555.1'] = 'repA'
cds_dict2['CAI2555830.1'] = '01250'
cds_dict2['CAI2561695.1'] = 'giant1' #02080
cds_dict2['CAI2561700.1'] = 'giant2' #02090
cds_dict2['CAI2561781.1'] = 'giant3' #02100
cds_dict2['CAI2561859.1'] = 'giant4' #02110
cds_dict2['CAI2561935.1'] = 'giant5' #02120

# =============================================================================
# 3. Process the dataframe for the first experiment
# =============================================================================

logging.info(f'Processing experiment data from {exp1}...')
for index, row in exp1_df.iterrows(): #Loop through one of the dfs

    logging.info(f'Add missing rows from {exp1} to {exp2}...')
    protein_id = row['Majority protein IDs'] #Get the protein IDs
    included = (exp2_df['Majority protein IDs'] == protein_id).any() #Check if the other df includes the protein IDs
    if not included: #If it is not in the other dataframe
        row_add = {'Protein IDs': row['Protein IDs'], #Create a row with protein IDs and fasta headers 
                   'Majority protein IDs': protein_id,
                   'Fasta headers': row['Fasta headers']}
        exp2_df = exp2_df._append(row_add, ignore_index = True) #Add the row to the other dataframe
        
    #Retrieve the names of columns with LFQs
    LFQ_cols = [column for column in exp1_df.columns if 'LFQ' in column]
        
    blank_val = row[f'LFQ intensity {len(LFQ_cols)}'] #Retrieve the value of the blank
    if blank_val != np.nan and blank_val > 0: #If the value of the blank is not missing or 0
        for ix in range(1, len(LFQ_cols)): #Loop through the LFQ columns
            exp1_df.loc[index, f'LFQ intensity {ix}'] -= blank_val #Subtract the blank value
        
    #Same as above, but for peptides
    pep_cols = [column for column in exp1_df.columns if 'Peptide' in column]
        
    blank_val = row[f'Peptides {len(pep_cols)}'] #Retrieve the value of the blank
    if blank_val != np.nan and blank_val > 0: #If the value of the blank is not missing or 0
        for ix in range(1, len(LFQ_cols)): #Loop through the LFQ columns
            exp1_df.loc[index, f'Peptides {ix}'] -= blank_val #Subtract the blank value
        
logging.info('Rename the columns to add sample names...')
columns = list(exp1_df.columns) #Retrieve dataframe columns
new_columns = [] #Create list to store new column names
for column in columns: #Loop through old columns
    column_end = column.split(' ')[-1] #Retrieve the last part of the column name
    column_rest = column.split(' ')[:-1] #Retrieve the rest of the column name
    new_column = column #Set the new column name to be the same as the old name
    
    #If the last part of the column name is a digit in the dictionary with sample names
    if column_end.isdigit() and int(column_end) in exp1_dict.keys():
        new_column_end = exp1_dict[int(column_end)] #Get the new end of the column name
        new_column = column_rest + [new_column_end] #Add it to the rest of the column information
        new_column = ' '.join(new_column) #Convert it back to a string
    #If the column name is numeric, but it is the blank
    elif column_end.isdigit(): 
        new_column_end = 'blank7' #Add "blank" to the end of the column name
        new_column = column_rest + [new_column_end]
        new_column = ' '.join(new_column)
        
    new_columns.append(new_column) #Add the new column name to a list
exp1_df.columns = new_columns #Replace the old column names with the new ones

#Sort the dataframe by protein ID and update the index
exp2_df.sort_values(['Majority protein IDs'], ignore_index = True, 
                    inplace = True)

# =============================================================================
# 4. Process the dataframe for the second experiment
# =============================================================================

logging.info(f'Processing experiment data from {exp2}...')
for index, row in exp2_df.iterrows(): #Same thing, but with the opposite dataframes
    logging.info(f'Add missing rows from {exp2} to {exp1}...')
    protein_id = row['Majority protein IDs']
    included = (exp1_df['Majority protein IDs'] == protein_id).any()
    if not included:
        row_add = {'Protein IDs': row['Protein IDs'], 
                   'Majority protein IDs': protein_id, 
                   'Fasta headers': row['Fasta headers']}
        exp1_df = exp1_df._append(row_add, ignore_index = True)

    #Retrieve the names of columns with LFQs
    LFQ_cols = [column for column in exp2_df.columns if 'LFQ' in column]
        
    blank_val = row[f'LFQ intensity {len(LFQ_cols)}'] #Retrieve the value of the blank
    if blank_val != np.nan and blank_val > 0: #If the value of the blank is not missing or 0
        for ix in range(1, len(LFQ_cols)): #Loop through the LFQ columns
            exp2_df.loc[index, f'LFQ intensity {ix}'] -= blank_val #Subtract the blank value
            
    #Same as above, but for peptides
    pep_cols = [column for column in exp2_df.columns if 'Peptide' in column][3:]
        
    blank_val = row[f'Peptides {len(pep_cols)}'] #Retrieve the value of the blank
    if blank_val != np.nan and blank_val > 0: #If the value of the blank is not missing or 0
        for ix in range(1, len(LFQ_cols)): #Loop through the LFQ columns
            exp2_df.loc[index, f'Peptides {ix}'] -= blank_val #Subtract the blank value
            
logging.info('Rename the columns to add sample names...')
columns = list(exp2_df.columns) #Retrieve dataframe columns
new_columns = [] #Create list to store new column names
for column in columns: #Loop through old columns
    column_end = column.split(' ')[-1] #Retrieve the last part of the column name
    column_rest = column.split(' ')[:-1] #Retrieve the rest of the column name
    new_column = column #Set the new column name to be the same as the old name
    
    #If the last part of the column name is a digit in the dictionary with sample names
    if column_end.isdigit() and int(column_end) in exp2_dict.keys():
        new_column_end = exp2_dict[int(column_end)] #Get the new end of the column name
        new_column = column_rest + [new_column_end] #Add it to the rest of the column information
        new_column = ' '.join(new_column) #Convert it back to a string
    #If the column name is numeric, but it is the blank
    elif column_end.isdigit(): 
        new_column_end = 'blank22' #Add "blank" to the end of the column name
        new_column = column_rest + [new_column_end]
        new_column = ' '.join(new_column)
        
    new_columns.append(new_column) #Add the new column name to a list
exp2_df.columns = new_columns #Replace the old column names with the new ones

#Sort the dataframe by protein ID and update the index
exp1_df.sort_values(['Majority protein IDs'], ignore_index = True, 
                    inplace = True)

# =============================================================================
# 5. Combine the results from the two experiments
# =============================================================================

logging.info('Merge the dataframes from the two experiments...')
exp_df = pd.merge(exp1_df, exp2_df, on = ['Majority protein IDs', 
                                          'Protein IDs', 'Fasta headers'], 
                  how = 'inner')

logging.info('Drop columns with blank data...')
exp_df.drop([col for col in exp_df.columns if 'blank' in col or 'counts' in col], 
            axis = 1, inplace = True)
        
logging.info('Add locus tags and annotations...')
exp_df['Locus tags'] = exp_df['Majority protein IDs'].map(prot2loctag)
exp_df['Original annotations'] = exp_df['Majority protein IDs'].map(cds_dict1)
exp_df['Refined annotations'] = exp_df['Majority protein IDs'].map(cds_dict2)

logging.info('Create a new dataframe with sorted columns...')
#Get the list of columns in the dataframe
cols = list(exp_df.columns)
#Reorder the columns
column_order = [cols[-3]] + cols[:2] + cols[-2:] + sorted(cols[2:-3])
#Create a new dataframe with reordered columns
df = exp_df[column_order]

logging.info('Write the dataframe to a tab file...')
df.to_csv(f'{outdir}/all_comparisons.tsv', header = df.columns, 
          index = None, sep = '\t', mode = 'w')

logging.info('Done!')
