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
from scipy.stats import ttest_ind
from Bio import SeqIO

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/proteomics'

log = f'{workdir}/logs/01alt-parse_proteomics.log'

with open(log, 'w') as logfile: #Overwrite log file
    logfile.write('')
            
# #Redirect stdout and stderr to log file
# sys.stdout = open(log, 'a')
# sys.stderr = open(log, 'a')

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

logging.info(f'Loading {exp2} dataset...')
exp2_df = pd.read_excel(exp2, sheet_name = 'All')
[exp2_df.insert(25+i, f'LFQ intensity {i}', np.nan) for i in range(1, 3)] #Insert missing samples as NaNs
exp2_df = exp2_df[exp2_df['Majority protein IDs'].notna()]
max_index = max(exp2_dict.keys()) + 1
exp2_df[f'LFQ intensity {max_index}'] = exp2_df[f'LFQ intensity {max_index}'].fillna(0) #Change the nan of the blank to 0s (otherwise, n - nan = nan)

# =============================================================================
# 3. Create a dictionary with annotations
# =============================================================================
print(annot1)
with open(annot1) as handle: #Open annotation file
    gbk = SeqIO.parse(handle, 'genbank') #Parse file as GenBank
    for record in gbk: #Loop through records in the file
        for cds in record.features: #Loop through features in each record
            if cds.type == 'CDS' and 'gene' in cds.qualifiers: #If the feature is a CDS and there is a gene name in the annotations
                cds_dict1[cds.qualifiers['protein_id'][0]] = cds.qualifiers['gene'][0] #Add the gene name to the dictionary, using the protein ID as a key
            elif cds.type == 'CDS' and cds.qualifiers['protein_id'][0] not in cds_dict1.keys(): #If the gene name is not in the annotations
                cds_dict1[cds.qualifiers['protein_id'][0]] = '-' #Use a minus symbol as gene name

print(cds_dict1)
                
with open(annot2) as handle:
    gbk = SeqIO.parse(handle, 'genbank')
    for record in gbk:
        for cds in record.features:
            if cds.type == 'CDS' and 'gene' in cds.qualifiers: #If the feature is a CDS and there is a gene name in the annotations
                cds_dict2[cds.qualifiers['locus_tag'][0]] = cds.qualifiers['gene'][0] #Add the gene name to the dictionary, using the protein ID as a key
            elif cds.type == 'CDS' and cds.qualifiers['locus_tag'][0] not in cds_dict2.keys(): #If the gene name is not in the annotations
                cds_dict2[cds.qualifiers['locus_tag'][0]] = '-' #Use a minus symbol as gene name
                
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
# 3. Make sure that all rows have the same index (add missing genes as NaNs)
# =============================================================================

print('Merging datasets...')
j = 0
for index, row in exp1_df.iterrows(): #Loop through one of the dfs
    protein_id = row['Majority protein IDs'] #Get the protein IDs
    included = (exp2_df['Majority protein IDs'] == protein_id).any() #Check if the other df includes the protein IDs
    if not included:
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

exp2_df.sort_values(['Majority protein IDs'], ignore_index = True, inplace = True) #Sort the dataframe by protein ID and update the index

for index, row in exp2_df.iterrows(): #Same thing, but with the opposite dataframes
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

exp1_df.sort_values(['Majority protein IDs'], ignore_index = True, 
                    inplace = True)

exp_df = pd.merge(exp1_df, exp2_df, on = ['Majority protein IDs', 
                                          'Protein IDs', 'Fasta headers'], 
                  how = 'inner')

exp_df.drop([col for col in exp_df.columns if ('LFQ' not in col or 'blank' in col) and ('ID' not in col)], axis = 1, 
            inplace = True)

# # =============================================================================
# # 4. Loop through the combined dataframe and subtract the blank for the samples
# # of each experiment
# # =============================================================================

# for index, row in exp_df.iterrows():
#     for ix in range(0, 7):
        

# # =============================================================================
# # 4. Create subsets using dictionaries and loops    
# # =============================================================================

# print('Creating subsets...')
# df_list = [] #Create an empty list to save dataframes

# i = 1 #Variable to loop
# while i < max(exp1_dict.keys()): #As long as i < the maximum index of the samples
#     print(f'Processing subset {exp1_dict[i][:-2]}')
#     df = exp1_df[['Majority protein IDs', 'Fasta headers', #Create a dataframe with the samples of interest only
#                   f'LFQ intensity {i}', f'LFQ intensity {i+1}', 
#                   f'LFQ intensity {i+2}']]
    
#     column_list = list(df.columns[:2]) + [f'LFQ_{exp1_dict[i]}', #Variable to update the name of the LFQ columns in the dataframe
#                                     f'LFQ_{exp1_dict[i+1]}', 
#                                     f'LFQ_{exp1_dict[i+2]}']
#     df.columns = column_list #Apply column name changes to dataframe
    
#     prot_ids = list(df['Majority protein IDs'])
#     labels = [cds_dict[pid] for pid in prot_ids]

#     df.insert(loc = 2, column = 'Gene names', value = labels)


#     df_list.append(df) #Add dataframe to the list
#     i += 3 #Increase the value of i by 3

# i = 1 #Same thing, but for the second dataset
# while i < max(exp2_dict.keys()):
#     print(f'Subset {exp2_dict[i][:-2]}')
#     df = exp2_df[['Majority protein IDs', 'Fasta headers',
#                   f'LFQ intensity {i}', f'LFQ intensity {i+1}', 
#                   f'LFQ intensity {i+2}']]
    
#     column_list = list(df.columns[:2]) + [f'LFQ_{exp2_dict[i]}', 
#                                     f'LFQ_{exp2_dict[i+1]}', 
#                                     f'LFQ_{exp2_dict[i+2]}']
    
#     df.columns = column_list
    
#     prot_ids = list(df['Majority protein IDs'])
#     labels = [cds_dict[pid] for pid in prot_ids]

#     df.insert(loc = 2, column = 'Gene names', value = labels)

    
#     df_list.append(df)
#     i += 3

# # =============================================================================
# # 5. Loop through subsets and calculate p-values when possible
# # =============================================================================

# print('Calculate p-values...')
# comparisons_list = []
# df_list.reverse() #Reverse the order of the dataframes
# for df1 in df_list: #Loop through all the dataframes with subsets of the data
#     cval_bool = df1.iloc[:, 3].name #Get the name of the sample (includes information about sample conditions)
#     print(f'Processing condition {cval_bool[:-2]}...')
#     cval = cval_bool.split('_') #Divide the information
#     c1 = cval[1][:3] #Condition 1 (mucoid - muc / inhibitor - inh)
#     c2 = cval[1][3] #Condition 2 (5% sucrose - S / 0.5% fructose - F)
#     c3 = cval[2][:-2] #Condition 3 (log phase - log / stationary phase - stat)
#     val_list = [c1, c2, c3] #Add conditions to a list
#     for df2 in df_list: #Loop through dataframes again
#         cval2 = df2.iloc[:, 3].name #Get the sample information from the other dataframe
#         no_sim = sum([val in cval2.replace('LFQ_', '') for val in val_list]) #Get the number of conditions that are similar between the dataframes
#         if cval_bool != cval2 and no_sim >= 2: #The dataframes have to differ, but at least two conditions have to be identical (dextranase disregarded)
#             print(f'Comparing {cval_bool} and {cval2}, with {no_sim} similar conditions')
#             comparison_df = pd.merge(df1, df2, on = ['Majority protein IDs', 
#                                                      'Fasta headers', 'Gene names']) #Merge the dataframes keeping the common columns
            
#             label1 = cval_bool[:-2].replace('LFQ_', '') #Create a label for one condition
#             label2 = cval2[:-2].replace('LFQ_', '') #Create a label for the other condition
#             labels = [label1, label2]
#             labels.sort(reverse = True)
#             if labels not in comparisons_list: #Check that the opposite comparison has not been done already
#                 comparisons_list.append(labels) #Add comparison to the dictionary
            
#                 dataset1 = list(comparison_df.iloc[:, 3:6].columns) #Get the data from the first condition
#                 dataset2 = list(comparison_df.iloc[:, 6:9].columns) #Get the data from the second condition
#                 print(dataset1, dataset2)
#                 comparison_df[f'avg_{label1}'] = comparison_df[dataset1].mean(axis=1) #Get the mean LFQ for one condition
#                 comparison_df[f'avg_{label2}'] = comparison_df[dataset2].mean(axis=1) #Get the mean LFQ for the other condition
#                 comparison_df.replace(0, np.nan, inplace = True) #Replace the averages of 0 with NaNs
#                 comparison_df[f'ratio_{label1}/{label2}'] = comparison_df[f'avg_{label1}']/comparison_df[f'avg_{label2}'] #Create a column with the ratio between the averages
#                 comparison_df[f'ratio_{label2}/{label1}'] = comparison_df[f'avg_{label2}']/comparison_df[f'avg_{label1}'] #Create a column with the opposite ratio

#                 # pval_list = [] #Create a list to store pvalues
#                 # for index, row in comparison_df.iterrows(): #Loop through the dataframe
#                 #     pval = ttest_ind(list(row[3:6].values), list(row[6:9].values), 
#                 #                      equal_var = False) #Calculate the p-value for each row
#                 #     pval_list.append(pval.pvalue) #Add the p-value to the list
 
#                 # comparison_df['p-value'] = pval_list #Create a column from the p-value list
                
#                 # comparison_df.dropna(axis = 0, thresh = 4, inplace = True) #Remove all the columns where all the LFQ values are NaNs
                
#                 # comparison_df.sort_values('p-value', ignore_index = True, #Sort the results by p-value
#                 #                           inplace = True)
                
#                 file_title = f'{label1}_vs_{label2}' #Set the title of the output file
#                 comparison_df.to_csv(f'{outdir}/{file_title}.tsv', #Write the dataframe to a tab file
#                                      header = comparison_df.columns, 
#                                      index = None, sep = '\t', mode = 'w')