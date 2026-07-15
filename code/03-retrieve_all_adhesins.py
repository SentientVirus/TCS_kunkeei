#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:11:46 2025

This script retrieves the genes that are equivalent to the adhesins from
H3B1-04J from the A. kunkeei genomes and creates a file with the locus tags
and another file with the sequences for each gene.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, pandas as pd
import re
from Bio import SeqIO
import logging, sys
import time

# =============================================================================
# 0. Logging
# =============================================================================

workdir = os.path.expanduser('~') + '/mucoid_project/adhesins' #Path to the working directory

logdir = f'{workdir}/logs' #Path to the log directory

#Create the log directory if it doesn't exist
if not os.path.exists(logdir):
    os.makedirs(logdir)
    
#Path to the log file
log = f'{logdir}/03-retrieve_all_adhesins.log'
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
# 1. Define input functions
# =============================================================================

def is_sublist(sublist, main_list):
    '''Function to check if a list is a sublist of another list.
    Inputs:
        sublist
        main_list'''
    sum_list = sum([element in main_list for element in sublist])
    if sum_list == len(sublist):
        return True
    else: return False

class adhesin:
    '''
    Class to store the locus tag, type and strain of each adhesin.
    '''
    
    def __init__(self, locus_tag, strain, adh_type = '?', ref_locus = '', 
                 pfams = [], domain_names = []): #Define class
        self.locus_tag = locus_tag #Locus tag
        self.ref_locus = ref_locus #Reference locus tag
        self.strain = strain #Strain
        self.type = adh_type #Adhesin type
        self.pfams = pfams #List of Pfam annotations
        self.domain_names = domain_names #List of domain names
        
    def get_type(self): #Function to assign a type to the adhesin based on Pfam annotations
        if self.pfams == ['PF05737']:
            self.type = 'collagen-binding'
            self.ref_locus = 'AKUH3B104J_00510'
        elif self.pfams == ['PF13632']:
            self.type = 'Gtf2'
            self.ref_locus = 'AKUH3B104J_00520'
        elif 'PF07564' in self.pfams or 'PF03382' in self.pfams: #and 'LDX55' not in self.locus_tag and (self.locus_tag.split('_')[1].startswith('01') or self.locus_tag.split('_')[1].startswith('RS006') or self.locus_tag in ['MUB42_02735', 'APS55_RS07870']):
            self.type = 'SH3b'
            loctag_no = self.locus_tag.split('_')[1] #Retrieve the number of the locus tag
            loctag_no = int(re.sub('[^0-9]', '', loctag_no)) #Convert it to integer
            if 'PF00746' in self.pfams: #Assign types based on Pfam domains if possible
                self.ref_locus = 'AKUH3B203M_02100'
            elif is_sublist(['PF03382', 'PF19087'], self.pfams) and loctag_no < 3000: #If not possible, use locus tag information too
                self.ref_locus = 'AKUH3B202X_01220'
            elif loctag_no < 10000:
                self.ref_locus = 'AKUH3B104J_01250'
            else: self.ref_locus = 'None'
        #Complex selection for adhesin 14310. Note that this adhesin is LPXTG 3-5 and can be further subdivided into types
        elif is_sublist(sorted(['PF00746', 'PF19258']), sorted(self.pfams)) or (is_sublist(sorted(['PF00746', 'PF17966']), sorted(self.pfams)) and 'PF17965' not in self.pfams and 'PF12799' not in self.pfams) or self.locus_tag == 'AKUH3B104J_14300':
            self.type = 'MubB2+LPXTG'
            self.ref_locus = 'AKUH3B104J_14310'
        #Selection for the adhesin LPXTG-8, which can be located in the plasmid (truncated in H3B1-04X)
        elif is_sublist(sorted(['PF17966', 'PF17965']), self.pfams) or self.locus_tag == 'AKUH3B104X_PLPX00300':
            self.type = 'MucBP+LPXTG'
            if ('PLPX' in self.locus_tag or int(self.locus_tag.replace('RS', '').split('_')[1]) < 1200 or 'MUB' in self.locus_tag) and self.locus_tag not in ['K2W83_RS00655', 'AAPFHON13_01040']:
                self.ref_locus = 'AKUH3B104J_01020'
            else:
                self.ref_locus = 'None'
    def __str__(self): #What print(class) returns
        return f'<Adhesin {self.locus_tag} from strain {self.strain} is {self.type}>'

# =============================================================================
# 2. Define input paramters
# =============================================================================

#Dictionary to color strains by phylogroup
phylogroup = {'A0901': 'C', 'A1001': 'F', 'A1002': 'C',
              'A1003': 'A', 'A1201': 'C', 'A1202': 'B',
              'A1401': 'B', 'A1404': 'E', 'A1802': 'C',
              'A1803': 'C', 'A1805': 'B', 'A2001': 'C',
              'A2002': 'F', 'A2003': 'F', 'A2101': 'B',
              'A2102': 'C', 'A2103': 'A', 'G0101': 'A', 
              'G0102': 'A', 'G0103': 'A', 'G0401': 'A',
              'G0402': 'A', 'G0403': 'B', 'G0404': 'A',
              'G0405': 'A', 'G0406': 'B', 'G0407': 'A',
              'G0408': 'A', 'G0410': 'A', 'G0412': 'A',
              'G0414': 'A', 'G0415': 'A', 'G0417': 'A',
              'G0420': 'B', 'G0601': 'A', 'G0602': 'A',
              'G0702': 'A', 'G0801': 'A', 'G0802': 'A',
              'G0803': 'A', 'G0804': 'A', 'Fhon2': 'A',
              'H1B1-04J': 'A', 'H1B1-05A': 'A', 'H1B3-02M': 'A',
              'H2B1-05J': 'C', 'H3B1-01A': 'A', 'H3B1-01J': 'A',
              'H3B1-01X': 'C', 'H3B1-02A': 'A', 'H3B1-02X': 'C',
              'H3B1-03J': 'A', 'H3B1-03X': 'C', 'H3B1-03M': 'A', 
              'H3B1-04J': 'A', 'H3B1-04X': 'A', 'H3B1-07A': 'A', 
              'H3B1-09M': 'A', 'H3B1-10M': 'A', 'H3B1-11A': 'A', 
              'H3B1-11M': 'A', 'H3B2-02M': 'C', 'H3B2-02X': 'A', 
              'H3B2-03J': 'A', 'H3B2-03M': 'C', 'H3B2-04J': 'A', 
              'H3B2-04M': 'C', 'H3B2-05J': 'A', 'H3B2-06M': 'C', 
              'H3B2-07X': 'A', 'H3B2-08X': 'A', 'H3B2-09X': 'B', 
              'H4B1-01A': 'C', 'H4B1-02A': 'C', 'H4B1-03J': 'C', 
              'H4B1-04A': 'C', 'H4B1-11J': 'C', 'H4B1-14J': 'C', 
              'H4B1-16J': 'C', 'H4B2-02J': 'A', 'H4B2-03M': 'C', 
              'H4B2-04J': 'A', 'H4B2-05J': 'A', 'H4B2-06J': 'C', 
              'H4B2-10M': 'C', 'H4B2-11M': 'A', 'H4B3-03J': 'C', 
              'H4B4-02J': 'A', 'H4B4-03J': 'C', 'H4B4-04J': 'C', 
              'H4B4-05J': 'C', 'H4B4-06M': 'C', 'H4B4-10M': 'C', 
              'H4B4-11M': 'C', 'H4B4-12M': 'A', 'H4B5-01J': 'A', 
              'H4B5-02X': 'A', 'H4B5-03X': 'A', 'H4B5-04J': 'B',
              'H4B5-05J': 'B', 'H4B5-07J': 'A', 'H4B5-07X': 'A', 
              'H4B5-08X': 'A', 'MP2': 'B', 'IBH001': 'C', 
              'DSMZ12361': 'A', 'HNS-8': 'X', 'GYUN-333': 'X', 
              'Fhon13': 'AAP'}

#Strains to exclude
exclude = ['G0403', 'G0406', 'G0420', 'H3B1-02X', 'H3B2-03M', 'H4B1-02A', 
           'H4B4-05J', 'H4B4-10M']

#Adhesin types
adhesins = ['MucBP+LPXTG', 'MubB2+LPXTG', 'Gtf2', 'collagen-binding', 'SH3b']

#Pfam IDs of the domains of interest
pfam_domains = ['PF06458', 'PF19087', 'PF05737', 'PF19258', #MucBP, DUF5776, collagen-binding, signal peptide
                'PF13632', 'PF17966', 'PF17965',  #gtf2, MucB2, MucBP_2
                'PF00746', 'PF07564', 'PF03382',  #LPXTG, EBH (DUF1542), DUF285
                'PF12799'] #Leucine-rich repeats (to exclude certain sequences)

#Pfam names of the domains of interest (abbreviated)
pfam_names = ['MucBP', 'DUF5776', 'collagen-binding', 'Gtf2', 'MucB2', 
              'MucBP_2', 'LPXTG']

#Representative strains
repr_strains = ['DSMZ12361', 'IBH001', 'GYUN-333', 'HNS-8', 'A0901', 
                'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805', 
                'Fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A', 
                'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 
                'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 'H4B2-06J', 'H4B4-02J', 
                'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 
                'H4B5-04J', 'H4B5-05J', 'MP2', 'Fhon13']

#Paths to inputs and outputs
workdir = os.path.expanduser('~') + '/mucoid_project/adhesins'
annot_dir = f'{workdir}/interproscan/locus_tags'
faa_dir = os.path.expanduser('~') + '/Akunkeei_files/faa'
prot2loctag = f'{workdir}/metadata/prot_id_loctag.tsv' #TSV with locus tags and protein IDs
outdir = f'{workdir}/results/adhesin_lists'
seqdir = f'{workdir}/sequences/adhesins'

#List of InterProScan files
interpro_files = [f'{annot_dir}/{file}' for file in os.listdir(annot_dir) if file.endswith('.tsv') and 'M-0' not in file]

#Create output directories if they don't exist
[os.makedirs(out) for out in [outdir, seqdir] if not os.path.exists(out)]

with open(prot2loctag) as handle: #Open TSV with locus tag and protein IDs
    #Create a dictionary mapping each locus tag to a protein ID
    loctag_dict = {line.split('\t')[1].strip(): line.split('\t')[0] for line in handle}
    
#Initialize variables
gene_list = [] 
curr_adh = adhesin('', '')
adh_dict = {}
strains = []
to_retrieve = {}
    
# =============================================================================
# 2. Loop through InterProScan annotations and retrieve adhesin information
# =============================================================================

logging.info('Parse InterProScan annotations to retrieve adhesins...')
for file in sorted(interpro_files): #Loop through InterProScan annotations
    with open(file) as ips: #Open file
        annot_df = pd.read_csv(ips, sep = '\t', header = None) #Parse annotations as a dataframe
        pfam_df = annot_df[annot_df[4].isin(pfam_domains)] #Filter out the lines that do not contain the domains of interest
        strain = os.path.basename(file).split('_')[0].replace('.tsv', '') #Retrieve the strain name
        strains.append(strain) #Add the strain to a list
        adh_dict[strain] = [] #Create an empty list of adhesins in the strain
        
        for index, row in pfam_df.iterrows(): #In the Pfam lines, search for any of the domains
            if row[0] not in gene_list: #If the row is not in the list of genes 
                #If any of the domains is in the protein, record this information in an object of class adhesin
                gene_list.append(row[0]) #Add the locus tag to the list of genes
                next_adh = adhesin(row[0], strain, pfams = [], domain_names = []) #Create adhesin object
                if next_adh.locus_tag != curr_adh.locus_tag: #If the next locus tag in the file in different from the current one
                    curr_adh.get_type() #Assign a type to the adhesin
                    if curr_adh.locus_tag != '': #If there is a current locus tag
                        adh_dict[curr_adh.strain].append(curr_adh) #Add it to the list of adhesins in the strain
                    logging.info(f'Adding {curr_adh.locus_tag} to the list of adhesins in {strain}!')
                    curr_adh = next_adh #Move on to the next adhesin
            curr_adh.pfams.append(row[4]) #Add the Pfam domains to the adhesin object
            curr_adh.domain_names.append(row[5]) #Add the Pfam domain names to the adhesin object
                
# =============================================================================
# 3. Loop through the adhesins types and write the sequences to a file        
# =============================================================================

for adhesin in adhesins: #Loop through the adhesin types
    outfile = f'{outdir}/{adhesin}_list.tsv' #Set path to the output file with a list of adhesins
    faa_outfile = f'{seqdir}/{adhesin}.faa' #Set the path to the output FASTA file
    faa_repset = f'{seqdir}/{adhesin}_repset.faa' #Set the path to the output FASTA file (representative strains only)
    
    logging.info(f'Processing adhesin {adhesin}...')
    logging.info('Output TSV: {outfile}\nOutput FASTA: {faa_outfile}, {faa_repset}')
    
    with open(outfile, 'w') as handle: #Open the out file
        if adhesin == 'MucBP+LPXTG': #Depending on the type, write different headers
            handle.write('Strain\tAdhesin_locus1\tAdhesin_locus2\n')
        elif adhesin == 'MubB2+LPXTG':
            handle.write('Strain\tAdhesin_locus1\tAdhesin_locus2\tAdhesin_locus3\n')
        elif adhesin == 'SH3b':
            handle.write('Strain\tSH3b_locus1\tLPXTG\tSH3b_locus2\tSH3b_locus3\n')
        else:
            handle.write('Strain\tAdhesin_locus\n')
            
        for strain in sorted(strains): #Loop through the strains
            #Select the adhesins in the strain that belong to the right type
            adh_type = [adh for adh in adh_dict[strain] if adh.type == adhesin]
            adh_type.sort(key=lambda x: x.locus_tag, reverse = False) #Sort the adhesins by locus tag
            to_retrieve[strain] = [adh.locus_tag for adh in adh_type] #Get a list of adhesin locus tags
            
            if adhesin == 'MucBP+LPXTG': #Retrieve the different adhesin subtypes
                adh1 = [adh for adh in adh_type if adh.ref_locus == 'AKUH3B104J_01020']
                adh2 = [adh for adh in adh_type if adh.ref_locus == 'None']
                
                if len(adh1) == 0: #If the first subtype is missing 
                    adh1_string = '-' #Indicate so
                else: #Otherwise
                    adh1_string = '' #Create an empty string
                    for adh in adh1: #Loop through the adhesin locus tags
                        adh1_string += adh.locus_tag + ', ' #Add the locus tags to the string, separated by commas
                    adh1_string = adh1_string[:-2] #Remove the comma and space at the end
                if len(adh2) == 0: #If the second subtype is missing
                    adh2_string = '-' #Indicate so
                else: #Otherwise, add the locus tags as before
                    adh2_string = ''
                    for adh in adh2:
                        adh2_string += adh.locus_tag + ', '
                    adh2_string = adh2_string[:-2]
                    
                #If at least one of the adhesins is found, write the locus tags to a file
                if not (adh1_string == '-' and adh2_string == '-'): 
                    handle.write(f'{strain}\t{adh1_string}\t{adh2_string}\n')
                    
            elif adhesin == 'MubB2+LPXTG': #If it is the LPXTG-3-5 adhesins
                rev = False #Set a variable indicating if the assembly of the strain is reversed
                if strain == 'MP2': #If the strain is MP2
                    rev = True #Set the variable to True
                    
                #Retrieve the list of adhesins and sort them by locus tag
                adhs = [adh for adh in adh_dict[strain] if adh.ref_locus == 'AKUH3B104J_14310']
                adhs.sort(key=lambda x: x.locus_tag, reverse = rev)

                if len(adhs) == 2: #Depending on the number of adhesins and strain, write the subtypes to a file
                    handle.write(f'{strain}\t{adhs[0].locus_tag}\t{adhs[1].locus_tag}\t-\n')
                elif len(adhs) > 0 and strain in exclude:
                    handle.write(f'{strain}\t{adhs[0].locus_tag}\t-\t-\n')
                elif len(adhs) > 0 and (phylogroup[strain] == 'F' or phylogroup[strain] == 'C'):
                    handle.write(f'{strain}\t-\t-\t{adhs[0].locus_tag}\n')
                elif len(adhs) > 0:
                    handle.write(f'{strain}\t-\t{adhs[0].locus_tag}\t-\n')
                    
            elif adhesin == 'SH3b': #Do the same for the SH3b as for the first adhesin, but with more subtypes
                adh1 = [adh for adh in adh_type if adh.ref_locus == 'AKUH3B104J_01250']
                adh2 = [adh for adh in adh_type if adh.ref_locus == 'AKUH3B202X_01220']
                adh3 = [adh for adh in adh_type if adh.ref_locus == 'None']
                LPXTG = [adh for adh in adh_type if adh.ref_locus == 'AKUH3B203M_02100']
                if len(adh1) == 0:
                    adh1_string = '-'
                else:
                    adh1_string = ''
                    for adh in adh1:
                        adh1_string += adh.locus_tag + ', '
                    adh1_string = adh1_string[:-2]
                if len(adh2) == 0:
                    adh2_string = '-'
                else:
                    adh2_string = ''
                    for adh in adh2:
                        adh2_string += adh.locus_tag + ', '
                    adh2_string = adh2_string[:-2]
                if len(adh3) == 0:
                    adh3_string = '-'
                else:
                    adh3_string = ''
                    for adh in adh3:
                        adh3_string += adh.locus_tag + ', '
                    adh3_string = adh3_string[:-2]
                if len(LPXTG) == 0:
                    LPXTG_string = '-'
                else:
                    LPXTG_string = ''
                    for adh in LPXTG:
                        LPXTG_string += adh.locus_tag + ', '
                    LPXTG_string = LPXTG_string[:-2]
                    
                handle.write(f'{strain}\t{adh1_string}\t{LPXTG_string}\t{adh2_string}\t{adh3_string}\n')
                
            else: #If the adhesin is not any of the previous types, just write the strain and locus tag
                [handle.write(f'{strain}\t{adh.locus_tag}\n') for adh in adh_type]
                
    logging.info(f'Writing {faa_outfile} and {faa_repset}!')
    #Open output files
    with open(faa_outfile, 'w') as faa_out, open(faa_repset, 'w') as faa_rep:
        for genome in sorted(list(to_retrieve.keys())): #Loop through strains
            faa = f'{faa_dir}/{genome}_protein.faa' #Retrieve the protein FASTA file
            check = False #Set the boolean to False
            #If it is one of the more recent strains
            if genome == 'HNS-8' or genome == 'GYUN-333':
                faa = faa.replace('/faa', '/new_genomes/faa') #Update the path
            if genome in repr_strains: #Loop through representative strains
                check = True #Set the boolean to True
            for record in SeqIO.parse(faa, 'fasta'): #Loop through records in the input FASTA
                if loctag_dict[record.id] in to_retrieve[genome]: #If the record is one of the adhesins to retrieve
                    record.id = loctag_dict[record.id] #Convert the protein ID to locus tag
                    record.description = '' #Remove the record description
                    SeqIO.write(record, faa_out, 'fasta') #Write the record to the general FASTA file
                    if check: #If the boolean is True (representative strain)
                        SeqIO.write(record, faa_rep, 'fasta') #Write the record to the FASTA file for representative strains

end_time = time.time() - start_time #Get total running time of the script
logging.info(f'Done! This script took {end_time/60:2f} minutes.')
            
        
