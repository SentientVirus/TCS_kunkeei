#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:23:32 2024

Script that reads all the Excel files with proteomics data, converts them into
tables and generates:
    
    1) A file with all the LFQ values of each sample
    2) One output file per comparison with the LFQs of both conditions,
       averages and p-values
       
The output files for each comparison should be used as inputs by other script
so that locus tags, signal peptide and TM helices predictions are added.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries
# =============================================================================

import numpy as np
import pandas as pd
import os
from scipy.stats import ttest_ind

# =============================================================================
# 1. Define paths to inputs and other variables
# =============================================================================

workdir = os.path.expanduser('~') + '/proteomics' #Working directory
outdir = f'{workdir}/files/parsed' #Output directory

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

exp1 = f'{workdir}/files/MS-24-030_MaxQuant_results.xlsx' #File with the results of the first MS dataset

exp2 = f'{workdir}/files/MS-24-038_MaxQuant_results.xlsx' #File with the results of the second MS dataset

exp1_dict = {1: 'mucSD_stat01', 2: 'mucSD_stat02', 3: 'mucSD_stat03', #Sample annotations of the first MS dataset
             4: 'inhS_stat04', 5: 'inhS_stat05', 6: 'inhS_stat06'}

exp2_dict = {1: 'mucS_log01', 2: 'mucS_log02', 3: 'mucS_log03', #Sample annotations of the second dataset
             4: 'mucF_log04', 5: 'mucF_log05', 6: 'mucF_log06',  
             7: 'mucF_stat04', 8: 'mucF_stat05', 9: 'mucF_stat06',
             10: 'mucFD_stat07', 11: 'mucFD_stat08', 12: 'mucFD_stat09',
             13: 'inhS_log10', 14: 'inhS_log11', 15: 'inhS_log12',
             16: 'inhF_log16', 17: 'inhF_log17', 18: 'inhF_log18',
             19: 'inhF_stat19', 20: 'inhF_stat20', 21: 'inhF_stat21'}

# =============================================================================
# Read input data
# =============================================================================

exp1_df = pd.read_excel(exp1, sheet_name = 'All') #Read results from MS-24-030 experiment
exp1_df = exp1_df.drop('i', axis = 1) #Remove column i
exp1_df = exp1_df[exp1_df['Majority protein IDs'].notna()] #Filter out empty rows and rows with column averages
max_index = max(exp1_dict.keys()) + 1 #Get the blank column index
exp1_df[f'LFQ intensity {max_index}'] = exp1_df[f'LFQ intensity {max_index}'].fillna(0) #Change the nan of the blank to 0s (otherwise, n - nan = nan)
for key in exp1_dict.keys(): #Loop through LFQ columns
    exp1_df[f'LFQ intensity {key}'] =  exp1_df[f'LFQ intensity {key}'] - exp1_df[f'LFQ intensity {max_index}'] #Substract the value of the blank from the column
    exp1_df[exp1_df[f'LFQ intensity {key}'] < 0] = np.nan #Convert values < blank to nan

exp2_df = pd.read_excel(exp2, sheet_name = 'All')
exp2_df['LFQ intensity 1'] = np.nan #Fill lost sample with NaNs
exp2_df['LFQ intensity 2'] = np.nan #Fill lost sample with NaNs
exp2_df = exp2_df[exp2_df['Majority protein IDs'].notna()]
max_index = max(exp2_dict.keys()) + 1
exp2_df[f'LFQ intensity {max_index}'] = exp2_df[f'LFQ intensity {max_index}'].fillna(0)
for key in exp2_dict.keys():
    exp2_df[f'LFQ intensity {key}'] =  exp2_df[f'LFQ intensity {key}'] - exp2_df[f'LFQ intensity {max_index}']

# =============================================================================
# Make sure that all rows have the same index (add missing genes as NaNs)
# =============================================================================

for index, row in exp1_df.iterrows(): #Loop through one of the dfs
    protein_id = row['Majority protein IDs'] #Get the protein IDs
    included = (exp2_df['Majority protein IDs'] == protein_id).any() #Check if the other df includes the protein IDs
    if not included: #If not
        row_add = {'Protein IDs': row['Protein IDs'], #Create a row with protein IDs and fasta headers 
                   'Majority protein IDs': protein_id, 
                   'Fasta headers': row['Fasta headers']}
        exp2_df = exp2_df._append(row_add, ignore_index = True) #Add the row to the other dataframe
        
exp2_df.sort_values(['Majority protein IDs'], ignore_index = True, inplace = True) #Sort the dataframe by protein ID and update the index

for index, row in exp2_df.iterrows(): #Same thing, but with the opposite dataframes
    protein_id = row['Majority protein IDs']
    included = (exp1_df['Majority protein IDs'] == protein_id).any()
    if not included:
        row_add = {'Protein IDs': row['Protein IDs'], 
                   'Majority protein IDs': protein_id, 
                   'Fasta headers': row['Fasta headers']}
        exp1_df = exp1_df._append(row_add, ignore_index = True)
        
exp1_df.sort_values(['Majority protein IDs'], ignore_index = True, inplace = True)

# =============================================================================
# Create subsets using dictionaries and loops    
# =============================================================================

df_list = [] #Create an empty list to save dataframes

i = 1 #Variable to loop
while i < max(exp1_dict.keys()): #As long as i < the maximum index of the samples
    df = exp1_df[['Protein IDs', 'Majority protein IDs', 'Fasta headers', #Create a dataframe with the samples of interest only
                  f'LFQ intensity {i}', f'LFQ intensity {i+1}', 
                  f'LFQ intensity {i+2}']]
    
    column_list = list(df.columns[:3]) + [f'LFQ_{exp1_dict[i]}', #Variable to update the name of the LFQ columns in the dataframe
                                    f'LFQ_{exp1_dict[i+1]}', 
                                    f'LFQ_{exp1_dict[i+2]}']
    df.columns = column_list #Apply column name changes to dataframe

    df_list.append(df) #Add dataframe to the list
    i += 3 #Increase the value of i by 3

i = 1 #Same thing, but for the second dataset
while i < max(exp2_dict.keys()):
    df = exp2_df[['Protein IDs', 'Majority protein IDs', 'Fasta headers', 
                  f'LFQ intensity {i}', f'LFQ intensity {i+1}', 
                  f'LFQ intensity {i+2}']]
    
    column_list = list(df.columns[:3]) + [f'LFQ_{exp2_dict[i]}', 
                                    f'LFQ_{exp2_dict[i+1]}', 
                                    f'LFQ_{exp2_dict[i+2]}']
    
    df.columns = column_list
    
    df_list.append(df)
    i += 3

# =============================================================================
# Loop through subsets and calculate p-values when possible
# =============================================================================

df_list.reverse() #Reverse the order of the dataframes
for df1 in df_list: #Loop through all the dataframes with subsets of the data
    cval_bool = df1.iloc[:, 3].name #Get the name of the sample (includes information about sample conditions)
    cval = cval_bool.split('_') #Divide the information
    c1 = cval[1][:3] #Condition 1 (mucoid - muc / inhibitor - inh)
    c2 = cval[1][3] #Condition 2 (5% sucrose - S / 0.5% fructose - F)
    c3 = cval[2][:-2] #Condition 3 (log phase - log / stationary phase - stat)
    val_list = [c1, c2, c3] #Add conditions to a list
    for df2 in df_list: #Loop through dataframes again
        cval2 = df2.iloc[:, 3].name #Get the sample information from the other dataframe
        no_sim = sum([val in cval2 for val in val_list]) #Get the number of conditions that are similar between the dataframes
        if cval_bool != cval2 and no_sim >= 2: #The dataframes have to differ, but at least two conditions have to be identical (dextranase disregarded)
            print(val_list, cval2, no_sim)
            comparison_df = pd.merge(df1, df2, on = ['Protein IDs', 
                                                     'Majority protein IDs', 
                                                     'Fasta headers']) #Merge the dataframes keeping the common columns
            
            dataset1 = list(comparison_df.iloc[:, 3:6].columns) #Get the data from the first condition
            dataset2 = list(comparison_df.iloc[:, 6:9].columns) #Get the data from the second condition
            label1 = cval_bool[:-2].replace('LFQ_', '') #Create a label for one condition
            label2 = cval2[:-2].replace('LFQ_', '') #Create a label for the other condition
            comparison_df[f'avg_{label1}'] = comparison_df[dataset1].mean(axis=1) #Get the mean LFQ for one condition
            comparison_df[f'avg_{label2}'] = comparison_df[dataset2].mean(axis=1) #Get the mean LFQ for the other condition
            comparison_df.replace(0, np.nan, inplace = True) #Replace the averages of 0 with NaNs
            comparison_df[f'ratio_{label1}/{label2}'] = comparison_df[f'avg_{label1}']/comparison_df[f'avg_{label2}'] #Create a column with the ratio between the averages

            pval_list = [] #Create a list to store pvalues
            for index, row in comparison_df.iterrows(): #Loop through the dataframe
                pval = ttest_ind(list(row[3:6].values), list(row[6:9].values)) #Calculate the p-value for each row
                pval_list.append(pval.pvalue) #Add the p-value to the list
            comparison_df['p-value'] = pval_list #Create a column from the p-value list
            
            comparison_df.dropna(axis = 0, thresh = 4, inplace = True) #Remove all the columns where all the LFQ values are NaNs
            
            comparison_df.sort_values('p-value', ignore_index = True, #Sort the results by p-value
                                      inplace = True)
            
            file_title = f'{label1}_vs_{label2}' #Set the title of the output file
            comparison_df.to_csv(f'{outdir}/{file_title}.tsv', #Write the dataframe to a tab file
                                 header = comparison_df.columns, 
                                 index = None, sep = '\t', mode = 'w')

    