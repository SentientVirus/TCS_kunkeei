#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 16:47:34 2025

Script to retrieve the sequences of the RR-TF and the histidine kinase.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
import time
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import logging, sys

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/adhesins' #Path to the working directory

logdir = f'{workdir}/logs' #Path to the log directory

#Create the log directory if it doesn't exist
if not os.path.exists(logdir):
    os.makedirs(logdir)
    
#Path to the log file
log = f'{logdir}/03a-get_RR-TF_HK_trees.log'
with open(log, 'w') as logfile: #Overwrite log file
    logfile.write('')
            
#Redirect stdout and stderr to log file
sys.stdout = open(log, 'a')
sys.stderr = open(log, 'a')

#Logging configuration
logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

start_time = time.time() #Get starting time

# =============================================================================
# 1. Define input variables
# =============================================================================

#Input directories
inpath = os.path.expanduser('~') + '/Akunkeei_files/gbff' #Path to the GenBank files
interpro_dir = f'{workdir}/interproscan/locus_tags' #Path to InterProScan annotations

#Output directories
outseqs = f'{workdir}/sequences/RR-TF_HK' #Path to the FASTA files with the sequences
outdir = outseqs.replace('sequences', 'alignments') #Path to the alignments
tree_dir = f'{workdir}/trees/RR-TF_HK' #Path to tree files

#Parameters/intermediate variables
threads = 8 #No. of threads to run software
gene_dict = {} #Dictionary to store gene information
#Representative strains to include in the analysis
repr_strains = ['DSMZ12361', 'IBH001', 'GYUN-333', 'HNS-8', 'A0901', 
                'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805', 
                'Fhon2', 'G0102', 'G0403', 'H1B1-04J', 'H1B1-05A', 
                'H1B3-02M', 'H3B1-11M', 'H3B1-04J', 'H3B1-04X', 'H3B1-03M', #The last one is an extra strain, added because of the plasmid gene
                'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 'H4B2-06J', 'H4B4-02J', 
                'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 
                'H4B5-04J', 'H4B5-05J', 'MP2', 'Fhon13']
new_strains = ['GYUN-333', 'HNS-8'] #Complete genomes that were published more recently

#Input files
#GenBank files of the more recently published genomes
extra_files = [inpath.replace('/gbff', '/new_genomes/gbff') + f'/{strain}_genomic.gbff' for strain in new_strains]
#Remaining GenBank files
infiles = [f'{inpath}/{file}' for file in sorted(os.listdir(inpath)) if file.endswith('.gbff') and 'M-0' not in file and any(rep in file for rep in repr_strains)]
infiles += extra_files #Merge the two lists
infiles = sorted(infiles) #Sort the input files

#Create output directories if they don't exist
new_dirs = [outseqs, outdir, tree_dir] #List of output directories
#List comprehension to create them if needed
[os.makedirs(new_dir) for new_dir in new_dirs if not os.path.exists(new_dir)]

# =============================================================================
# 2. Retrieve the RR-TF and HK sequences
# =============================================================================

for file in infiles: #Loop through input files
    strain = os.path.basename(file).replace('_genomic.gbff', '') #Retrieve strain name
    interpro_file = f'{interpro_dir}/{strain}.tsv' #Set the path to the InterProScan file
    
    with open(interpro_file) as handle: #Open InterProScan file
        annot = pd.read_csv(handle, sep = '\t', header = None) #Read annotations as dataframe
        #Set dataframe column names
        annot.columns = ['locus_tag', 'MD5_digest', 'length', 'analysis',
                         'analysis_accession', 'analysis_description', 
                         'start', 'end', 'score', 'status', 'date', 
                         'InterPro_accession', 'InterPro_description']
    #Retrieve the RR by annotation
    RR = annot[annot['analysis_description'].str.contains('ResD')]['locus_tag']
    if len(RR) > 0: #If the RR is present
        RR = RR.to_string().split('    ')[1] #Retrieve the locus tag
        logging.info(f'Retrieving RR-TF for strain {strain}...')
    else: 
        RR = None #Otherwise, set it to None
        logging.info(f'RR-TF not found in strain {strain}!')
    
    #Do the same for the histidine kinase
    HK = annot[annot['analysis_accession'] == 'G3DSA:3.30.565.10:FF:000013']['locus_tag']
    if len(HK) > 0: #If the HK is found
        HK = HK.to_string().split('    ')[1] #Retrieve the locus tag as a string
                                          
    #If the strain comes from the dataset from Dyrhage et al. (2022) and the RR is found
    elif strain not in ['DSMZ12361', 'HNS-8', 'GYUN-333', 'IBH001', 'MP2'] and RR != None:
        tag = int(RR.split('_')[1]) + 10 #Retrieve the HK as the locus tag after the RR
        HK = RR.split('_')[0] + '_0' + str(tag) #Convert back to string
        logging.info(f'Retrieving HK for strain {strain}...')
    #For certain other strains
    elif strain in ['HNS-8', 'GYUN-333', 'IBH001'] and RR != None:
        tag = int(RR.split('_')[1]) + 5 #Retrieve the HK as the locus tag after the RR
        HK = RR.split('_')[0] + '_' + str(tag) #Convert back to string
        logging.info(f'Retrieving HK for strain {strain}...')
    elif strain == 'DSMZ12361': #For the type strain
        #Retrieve the HK from a different previous annotation
        HK = annot[annot['analysis_accession'] == 'G3DSA:3.20.20.70:FF:000424']['locus_tag'].to_string().split('    ')[1]
        HK = HK.split('_')[0] + '_RS0' + str(int(HK.split('_RS')[1]) + 10) #Get the next locus tag
        logging.info(f'Retrieving HK for strain {strain}...')
    elif strain == 'MP2': #If the strain is MP2
        tag = int(RR.split('_')[1]) - 5 #Retrieve the HK as the locus tag before the RR (assembly on the reverse strand)
        HK = RR.split('_')[0] + '_' + str(tag) #Convert back to string
        logging.info(f'Retrieving HK for strain {strain}...')
    else: 
        HK = None #Otherwise, set the HK to None
        logging.info(f'HK not found in strain {strain}!')
    
    #If the RR or the HK are retrieved, assign gene types to each locus tag in the dictionary
    if RR != None:
        gene_dict[RR] = 'RR'
    if HK != None:
        gene_dict[HK] = 'HK'
        
    logging.info(f'{strain}, RR-TF {RR}, HK: {HK}') #Print retrieved locus tags
        
    with open(file) as gbff: #Open the GenBank file
        gbk = SeqIO.parse(gbff, 'genbank') #Parse the GenBank file
        
        for record in gbk: #Loop through records in the file
            for feature in record.features: #Loop through features in the records
                #If the feature is a CDS with an assigned locus tga
                if feature.type == 'CDS' and 'locus_tag' in feature.qualifiers.keys():
                    loctag = feature.qualifiers['locus_tag'][0] #Retrieve the locus tag
                    if loctag in gene_dict.keys(): #If the locus tag is in the dictionary with gene types
                        seq = feature.qualifiers['translation'][0] #Retrieve the gene sequence
                        #Create a new record where the ID is the locus tag and the name is the gene type
                        new_record = SeqRecord(Seq(seq), id = loctag, 
                                                name = gene_dict[loctag],
                                                description = '')
                        if file == infiles[0]: #If it is the first file
                            #Create/overwirte the output file
                            with open(f'{outseqs}/{gene_dict[loctag]}.faa', 'w') as rep_faa:
                                SeqIO.write(new_record, rep_faa, 'fasta') #Write the record to file
                        else: #If it is not the first file
                            #Append the next record to the file
                            with open(f'{outseqs}/{gene_dict[loctag]}.faa', 'a') as rep_faa:
                                SeqIO.write(new_record, rep_faa, 'fasta')
                                
# =============================================================================
# 3. Align the sequences wth MAFFT and generate a phylogeny with IQtree
# =============================================================================
                            
#Use the unique values in the dictionary (RR and HK) to loop through the output FASTA files
for file in [f'{outseqs}/{gene}.faa' for gene in list(set(gene_dict.values()))]:
    logging.info(f'Running MAFFT and IQtree on {file}...')
    #Define the path to the sequence alignment
    outfile = file.replace('.faa', '.mafft.faa').replace('sequences', 'alignments')
    #Run the sequence alignment with MAFFT
    subprocess.run(f'mafft-linsi --thread {threads} {file} > {outfile} 2>> {log};',
                    shell = True)
    #Create a phylogeny with IQtree
    subprocess.run(f'iqtree -nt AUTO -ntmax {threads} -redo -s {outfile} -st AA -msub nuclear -bb 1000 -bnni >> {log}', 
                    shell = True)
    #Move the IQtree output to the desired directory
    subprocess.run(f'mv {outfile}.* {tree_dir}', shell = True)
    
end_time = time.time() - start_time #Get total running time of the script
logging.info(f'This script took {end_time/60:2f} minutes.')
