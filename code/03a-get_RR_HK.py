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
import re
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
# 0. Function to run EMBOSS Needle
# =============================================================================

def needle_align_code(query_seq, target_seq):
    '''Function to perform a Needle alignment between a pair of sequences
    and retrieve the percentage of identity
    
    Inputs:
        - query_seq (str): One of the two sequences to compare
        - target_seq (str): The second sequence'''
        
    aseq = 'asis:' + query_seq #First sequence (must be processed as-is)
    bseq = 'asis:' + target_seq #Second sequence (also read as-is)
        
    #Define the command to run
    command = f'needle -asequence {aseq} -bsequence {bseq}'
    command += ' -gapopen 10 -gapextend 0.5 -aformat simple -outfile stdout'
    command += f' -snucleotide1 -snucleotide2 2>> {log}'
    
    #Run the command and save stdout and sterr to variables
    out_data = subprocess.run(command, shell = True, stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE, text = True)
    
    with open(log, 'a') as handle: #Open log file
        handle.write(out_data.stderr) #Write stderr to log file
        
    return out_data #Return the captured output

# =============================================================================
# 1. Define input variables
# =============================================================================

#Input directories
inpath = os.path.expanduser('~') + '/Akunkeei_files/cds' #Path to the CDS files
interpro_dir = f'{workdir}/interproscan/locus_tags' #Path to InterProScan annotations

#Output directories
outseqs = f'{workdir}/sequences/RR-TF_HK' #Path to the FASTA files with the sequences
outdir = outseqs.replace('sequences', 'alignments') #Path to the alignments
needle_outdir = f'{outdir}/needle' #Path to Needle results

#Parameters/intermediate variables
threads = 8 #No. of threads to run software
gene_dict = {} #Dictionary to store gene information
perc_ident = {} #Dictionary to store gene % of ID
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
#CDS files of the more recently published genomes
extra_files = [inpath.replace('/cds', '/new_genomes/fna') + f'/{strain}_cds_from_genomic.fna' for strain in new_strains]
#Remaining CDS files
infiles = [f'{inpath}/{file}' for file in sorted(os.listdir(inpath)) if file.endswith('.fna') and 'M-0' not in file and any(rep in file for rep in repr_strains)]
infiles += extra_files #Merge the two lists
infiles = sorted(infiles) #Sort the input files

#Retrieve also all the CDS files
all_infiles =  [f'{inpath}/{file}' for file in sorted(os.listdir(inpath)) if file.endswith('.fna') and 'M-0' not in file] + extra_files

#Create output directories if they don't exist
new_dirs = [outseqs, outdir, needle_outdir] #List of output directories
#List comprehension to create them if needed
[os.makedirs(new_dir) for new_dir in new_dirs if not os.path.exists(new_dir)]

p = re.compile("\((.*)\)") #Regular expression to parse Needle results

# =============================================================================
# 2. Retrieve the RR-TF and HK sequences
# =============================================================================

for file in all_infiles: #Loop through input files
    strain = os.path.basename(file).replace('_cds_from_genomic.fna', '') #Retrieve strain name
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
        
    with open(file) as cds: #Open the GenBank file
        cds_seqs = SeqIO.parse(cds, 'fasta') #Parse the GenBank file
        
        for record in cds_seqs: #Loop through records in the file
            #If the feature is a CDS with an assigned locus tga
            if 'locus_tag' in record.description:
                loctag = [info for info in record.description.split(']') if 'locus_tag=' in info][0]
                loctag = loctag.split('[locus_tag=')[-1]
                if loctag in gene_dict.keys(): #If the locus tag is in the dictionary with gene types
                    seq = record.seq #Retrieve the gene sequence
                    #Create a new record where the ID is the locus tag and the name is the gene type
                    new_record = SeqRecord(Seq(seq), id = loctag, 
                                            name = gene_dict[loctag],
                                            description = '')
                    if file == infiles[0]: #If it is the first file (representative strains)
                        #Create/overwrite the output file with representative sequences
                        with open(f'{outseqs}/{gene_dict[loctag]}_repset.fna', 'w') as rep_fna:
                            SeqIO.write(new_record, rep_fna, 'fasta') #Write the record to file
                    if file == all_infiles[0]: #If it is the first file (all strains)
                        #Create/overwirte the output file with all the sequences
                        with open(f'{outseqs}/{gene_dict[loctag]}_all.fna', 'w') as all_fna:
                            SeqIO.write(new_record, all_fna, 'fasta') #Write the record to file
                    if file != infiles[0] and file in infiles: #If it is not the first file (representative strains)
                        #Append the next record to the file
                        with open(f'{outseqs}/{gene_dict[loctag]}_repset.fna', 'a') as rep_fna:
                            SeqIO.write(new_record, rep_fna, 'fasta')
                    if file != infiles[0]: #If is not the first file (all strains)
                        #Append the next record to the file
                        with open(f'{outseqs}/{gene_dict[loctag]}_all.fna', 'a') as all_fna:
                            SeqIO.write(new_record, all_fna, 'fasta')
                                
# =============================================================================
# 3. Align the sequences wth MAFFT and generate a phylogeny with IQtree
# =============================================================================
                            
#Use the unique values in the dictionary (RR and HK) to loop through the output FASTA files
for file in [f'{outseqs}/{gene}_{suffix}.fna' for gene in list(set(gene_dict.values())) for suffix in ['repset', 'all']]:
    logging.info(f'Running MAFFT and IQtree on {file}...')
    #Define the path to the sequence alignment
    outfile = file.replace('.fna', '.mafft.fna').replace('sequences', 'alignments')
    
    #Define the path to summarized Needle results
    tab_outfile = outfile.replace('.mafft.fna', '_needle.tsv')
    with open(tab_outfile, 'w') as handle: #Open output file (create/overwrite mode)
        handle.write('Locus1\tLocus2\tPercent_ID\n') #Write headers
        
    #Run the sequence alignment with MAFFT
    subprocess.run(f'mafft-linsi --thread {threads} {file} > {outfile} 2>> {log};',
                    shell = True)
    
    records = [record for record in SeqIO.parse(file, 'fasta')]
    for i in range(0, len(records)-1):
        id1 = records[i].id
        seq1 = records[i].seq
        for j in range(i+1, len(records)):
            id2 = records[j].id
            seq2 = records[j].seq
            ident_out = needle_align_code(seq1, seq2)
            out_split = ident_out.stdout.split('\n') #Divide the output into lines
            ident = p.search(out_split[26]).group(1).replace('%', '') #Retrieve the percentage of identity
            perc_ident[(id1, id2)] = ident
            
            needle_out = f'{needle_outdir}/{id1}_vs_{id2}.txt'
            
            with open(needle_out, 'w') as handle:
                handle.write(ident_out.stdout)
                
            with open(tab_outfile, 'w') as handle:
                handle.write(f'{id1}\t{id2}\t{perc_ident}\n')
    
    # #Create a phylogeny with IQtree
    # subprocess.run(f'iqtree -nt AUTO -ntmax {threads} -redo -s {outfile} -st AA -msub nuclear -bb 1000 -bnni >> {log}', 
    #                 shell = True)
    # #Move the IQtree output to the desired directory
    # subprocess.run(f'mv {outfile}.* {tree_dir}', shell = True)
    
end_time = time.time() - start_time #Get total running time of the script
logging.info(f'This script took {end_time/60:2f} minutes.')
