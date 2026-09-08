#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:23:32 2024

Script that reads all the Excel files with proteomics data, converts them into
tables and generates:
    
    1) A file with all the LFQ values of each sample
    2) One output file per comparison with the LFQs of both conditions, and
    averages
       
The output files for each comparison should be used as inputs by other script
so that signal peptide and TM helices predictions are added.

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
# 1. Define paths to inputs and other variables
# =============================================================================

outdir = snakemake.params.outdir #Output directory

general_outfile = snakemake.output.general
print(general_outfile)

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

exp1 = snakemake.input.exp1 #File with the results of the first MS dataset

exp2 = snakemake.input.exp2 #File with the results of the second MS dataset

annot1 = snakemake.input.gbk1 #Reference GenBank file (original annotations)
annot2 = snakemake.input.gbk2 #Reference GenBank file (refined annotations)

exp1_dict = {1: 'mucSD_stat01', 2: 'mucSD_stat02', 3: 'mucSD_stat03', #Sample annotations of the first MS dataset
             4: 'aggS_stat04', 5: 'aggS_stat05', 6: 'aggS_stat06'}

exp2_dict = {1: 'mucS_log01', 2: 'mucS_log02', 3: 'mucS_log03', #Sample annotations of the second dataset
             4: 'mucF_log04', 5: 'mucF_log05', 6: 'mucF_log06',  
             7: 'mucF_stat07', 8: 'mucF_stat08', 9: 'mucF_stat09',
             10: 'mucFD_stat10', 11: 'mucFD_stat11', 12: 'mucFD_stat12',
             13: 'aggS_log13', 14: 'aggS_log14', 15: 'aggS_log15',
             16: 'aggF_log16', 17: 'aggF_log17', 18: 'aggF_log18',
             19: 'aggF_stat19', 20: 'aggF_stat20', 21: 'aggF_stat21'}

cds_dict1 = {} #Create an empty dictionary to store the original annotations
cds_dict2 = {} #Same, but for the modified annotations
prot2loctag = {} #Dictionary to assign locus tags to the protein IDs

# =============================================================================
# 2. Read input data
# =============================================================================

#Get experiment names
exp1_name = os.path.basename(exp1.replace('_results.xlsx', ''))
exp2_name = os.path.basename(exp2.replace('_results.xlsx', ''))

logging.info(f'Loading {exp1_name} dataset...')
exp1_df = pd.read_excel(exp1, sheet_name = 'All') #Read results from MS-24-030 experiment
exp1_df = exp1_df.drop('i', axis = 1) #Remove column i
exp1_df = exp1_df[exp1_df['Majority protein IDs'].notna()] #Filter out empty rows and rows with column averages
max_index = max(exp1_dict.keys()) + 1 #Get the blank column index
exp1_df[f'LFQ intensity {max_index}'] = exp1_df[f'LFQ intensity {max_index}'].fillna(0) #Change the nan of the blank to 0s (otherwise, n - nan = nan)
exp1_df[f'Peptides {max_index}'] = exp1_df[f'Peptides {max_index}'].fillna(0)

logging.info(f'Loading {exp2_name} dataset...')
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

logging.info(f'Processing experiment data from {exp1_name}...')
logging.info(f'Subtract blank values from {exp1_name}...\nAdd missing rows from {exp1_name} to {exp2_name}...')
for index, row in exp1_df.iterrows(): #Loop through one of the dfs

    protein_id = row['Majority protein IDs'] #Get the protein IDs
    included = (exp2_df['Majority protein IDs'] == protein_id).any() #Check if the other df includes the protein IDs
    if not included: #If it is not in the other dataframe
        row_add = {'Majority protein IDs': protein_id,
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
#Remove the Protein IDs column
exp1_df.drop('Protein IDs', axis = 1, inplace = True)

#Sort the dataframe by protein ID and update the index
exp2_df.sort_values(['Majority protein IDs'], ignore_index = True, 
                    inplace = True)

# =============================================================================
# 4. Process the dataframe for the second experiment
# =============================================================================

logging.info(f'Processing experiment data from {exp2_name}...')
logging.info(f'Subtract blank values from {exp2_name}...\nAdd missing rows from {exp1_name} to {exp2_name}...')
for index, row in exp2_df.iterrows(): #Same thing, but with the opposite dataframes
    protein_id = row['Majority protein IDs']
    included = (exp1_df['Majority protein IDs'] == protein_id).any()
    if not included:
        row_add = {'Majority protein IDs': protein_id, 
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
#Remove the Protein IDs column
exp2_df.drop('Protein IDs', axis = 1, inplace = True)

#Sort the dataframe by protein ID and update the index
exp1_df.sort_values(['Majority protein IDs'], ignore_index = True, 
                    inplace = True)

# =============================================================================
# 5. Combine the results from the two experiments
# =============================================================================

logging.info('Merge the dataframes from the two experiments...')
exp_df = pd.merge(exp1_df, exp2_df, on = ['Majority protein IDs', 
                                          'Fasta headers'], 
                  how = 'inner')

logging.info('Drop columns with blank data...')
exp_df.drop([col for col in exp_df.columns if 'blank' in col or 'counts' in col], 
            axis = 1, inplace = True)
        
logging.info('Add locus tags and annotations...')
exp_df['Majority locus tags'] = exp_df['Majority protein IDs'].map(prot2loctag)
exp_df['Original annotations'] = exp_df['Majority protein IDs'].map(cds_dict1)
exp_df['Refined annotations'] = exp_df['Majority protein IDs'].map(cds_dict2)

logging.info('Create a new dataframe with sorted columns...')
#Get the list of columns in the dataframe
cols = list(exp_df.columns)
#Reorder the columns
column_order = [cols[-3]] + cols[:2] + cols[-2:] + sorted(cols[3:-3])

#Create a new dataframe with reordered columns
df = exp_df[column_order]

logging.info('Write the dataframe to a tab file...')
df.to_csv(general_outfile, header = df.columns, 
          index = None, sep = '\t', mode = 'w')

logging.info('Done!')

# =============================================================================
# 6. Create subsets using dictionaries and loops    
# =============================================================================

logging.info('Creating subsets...')
df_list = [] #Create an empty list to save dataframes
column_list = list(df.columns)
info_columns = column_list[:5]
data_columns = [column for column in column_list if 'LFQ' in column]
print(info_columns, data_columns)

i = 0 #Variable to loop
while i < len(data_columns): #As long as i < the maximum index of the samples

    condition = ''.join([cx for cx in data_columns[i].split(' ')[-1] if not cx.isdigit()])
    logging.info(f'Processing subset {condition}')
    
    columns_select = data_columns[i:i+3]

    dfn = df[info_columns + columns_select]

    df_list.append(dfn) #Add dataframe to the list
    i += 3 #Increase the value of i by 3

# =============================================================================
# 5. Loop through subsets and calculate p-values when possible
# =============================================================================

logging.info('Calculate p-values...')
comparisons_list = []
df_list.reverse() #Reverse the order of the dataframes
for j in range(len(df_list)-1): #Loop through all the dataframes with subsets of the data
    df1 = df_list[j]
    cval_bool = df1.iloc[:, 7].name.split(' ')[-1] #Get the name of the sample (includes information about sample conditions)
    logging.info(f'Processing condition {cval_bool[:-2]}...')
    cval = cval_bool.split('_') #Divide the information
    c1 = cval[0][:3] #Condition 1 (mucoid - muc / aggregating - agg)
    c2 = cval[0][3] #Condition 2 (5% sucrose - S / 0.5% fructose - F)
    c3 = cval[1][:-2] #Condition 3 (log phase - log / stationary phase - stat)
    val_list = [c1, c2, c3] #Add conditions to a list
    
    logging.info(f'Phenotype: {c1}\tMedium: {c2}\tGrowth phase: {c3}')

    for k in range(j+1, len(df_list)): #Loop through dataframes again
        df2 = df_list[k]
        cval2_bool = df2.iloc[:, 7].name.split(' ')[-1] #Get the sample information from the other dataframe
        logging.info(f'Processing second condition {cval2_bool[:-2]}...')
        cval2 = cval2_bool.split('_')
        c4 = cval2[0][:3] #Condition 1 (mucoid - muc / aggregating - agg)
        c5 = cval2[0][3] #Condition 2 (5% sucrose - S / 0.5% fructose - F)
        c6 = cval2[1][:-2] #Condition 3 (log phase - log / stationary phase - stat)
        val2_list = [c4, c5, c6] #Add conditions to a list
        
        logging.info(f'Phenotype: {c4}\tMedium: {c5}\tGrowth phase: {c6}')
        
        no_sim = sum([val in val2_list for val in val_list]) #Get the number of conditions that are similar between the dataframes
        if cval_bool != cval2 and no_sim >= 2: #The dataframes have to differ, but at least two conditions have to be identical (dextranase disregarded)
            logging.info(f'Comparing {c1}{c2}_{c3} and {c4}{c5}_{c6}, with {no_sim} similar conditions')
            comparison_df = pd.merge(df1, df2, on = ['Majority locus tags',
                                                     'Majority protein IDs',
                                                     'Fasta headers', 
                                                     'Original annotations',
                                                     'Refined annotations']) #Merge the dataframes keeping the common columns
            
            label1 = cval_bool[:-2].replace('LFQ intensity ', '') #Create a label for one condition
            label2 = cval2_bool[:-2].replace('LFQ intensity ', '') #Create a label for the other condition
            labels = [label1, label2]
            labels.sort()
            if labels not in comparisons_list: #Check that the opposite comparison has not been done already
                comparisons_list.append(labels) #Add comparison to the dictionary
            
                dataset1 = list(comparison_df.iloc[:, 5:8].columns) #Get the data from the first condition
                dataset2 = list(comparison_df.iloc[:, 8:11].columns) #Get the data from the second condition
                logging.info(f'Dataset 1: {dataset1}\nDataset 2: {dataset2}')
                comparison_df.replace(0, np.nan, inplace = True) #Replace the averages of 0 with NaNs
                comparison_df[f'avg_{label1}'] = comparison_df[dataset1].mean(axis=1) #Get the mean LFQ for one condition
                comparison_df[f'avg_{label2}'] = comparison_df[dataset2].mean(axis=1) #Get the mean LFQ for the other condition
                comparison_df[f'ratio_{label1}/{label2}'] = comparison_df[f'avg_{label1}']/comparison_df[f'avg_{label2}'] #Create a column with the ratio between the averages
                comparison_df[f'ratio_{label2}/{label1}'] = comparison_df[f'avg_{label2}']/comparison_df[f'avg_{label1}'] #Create a column with the opposite ratio

                
                file_title = f'{labels[0]}_vs_{labels[1]}' #Set the title of the output file
                comparison_df.to_csv(f'{outdir}/{file_title}.tsv', #Write the dataframe to a tab file
                                     header = comparison_df.columns, 
                                     index = None, sep = '\t', mode = 'w')

    
