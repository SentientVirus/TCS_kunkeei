#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 17:11:34 2025

Script to retrieve the intergenetic regions between the MucBP adhesins, the
TF and other genes to search for a promoter.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import as_fasta
import subprocess
from pymsaviz import MsaViz

# =============================================================================
# 1. Define function to retrieve strain names
# =============================================================================

def get_unique_strains(tag_list):
    strains = []
    for tag in tag_list:
        strain = tag.split('_')[0].replace('AKU', '')
        if strain.startswith('H'):
            strain = strain[:4] + '-' + strain[4:]
        strain = strain.replace('K2W83', 'DSMZ12361').replace('MUB42', 'HNS-8').replace('FHON', 'Fhon').replace('APS55', 'MP2')
        if strain not in strains:
            strains.append(strain)
    strains = sorted(strains)
    return strains

# =============================================================================
# 2. Define global inputs
# =============================================================================

indir = os.path.expanduser('~') + '/Akunkeei_files/gbff' #Path to GenBanks
indir2 = indir.replace('/gbff', '/new_genomes/gbff') #Path to GenBanks of newest genomes
outdir = os.path.expanduser('~') + '/mucoid_project/adhesins/sequences/intergenic' #Output directory for FASTA files
plotdir = os.path.expanduser('~') + '/mucoid_project/adhesins/plots/intergenic' #Output directory for plots
consensus_file = f'{outdir}/consensus.fna' #Output consensus sequences
mafft_consensus = consensus_file.replace('.fna', '.mafft.fna') #Output consensus alignment
consensus_plot = f'{plotdir}/consensus.png' #Output consensus alignment plot

#Create output directories if they don't exist
[os.makedirs(dirn) for dirn in [outdir, plotdir] if not os.path.exists(dirn)]

with open(consensus_file, 'w') as handle: #Open the consensus file in write mode
    handle.write('') #Create it/Overwrite it

#Note: For the MubB2+LPXTG and the Gtf2&collagen-binding, that are placed next to each other in the genome, only the region upstream of the first locus is retrieved
adhesins = ['MucBP+LPXTG', 'MubB2+LPXTG', 'Gtf2'] #SH3b can also be added?

for adhesin in adhesins: #Loop through adhesins

# =============================================================================
# 2. Define the locus tags for the genes of interest (before/after the gene)
# =============================================================================

#This is for the regulator just 
# pos_tetR_tags = ['K2W83_RS00580', 'AKUFHON2_01070', #'AAPFHON13_00970', 
#               'AKUG0101_01080', 'AKUG0102_01070', 'AKUG0103_01070', 
#               'AKUG0401_01070', 'AKUG0402_01070', 'AKUG0403_PLPX00270',
#               'AKUG0404_01070', 'AKUG0405_01070', 'AKUG0406_PLPX00280',
#               'AKUG0407_01070', 'AKUG0408_01070', 'AKUG0410_01110',
#               'AKUG0412_01110', 'AKUG0414_01070', 'AKUG0415_01070',
#               'AKUG0417_01100', 'AKUG0420_PLPX00310', 'AKUG0601_01070',
#               'AKUG0602_01070', 'AKUG0702_01070', 'AKUG0801_01070',
#               'AKUG0802_01070', 'AKUG0803_01070', 'AKUG0804_01070',
#               'AKUH1B104J_01070', 'AKUH1B105A_00980', 'AKUH3B101A_01060',
#               'AKUH3B101J_01040', 'AKUH3B102A_01060', 'AKUH3B103J_01060',
#               'AKUH3B103M_PLPX00270', 'AKUH3B104J_01040', 'AKUH3B104X_PLPX00270',
#               'AKUH3B107A_01060', 'AKUH3B109M_01060', 'AKUH3B110M_01060',
#               'AKUH3B111A_PLPX00260', 'AKUH3B111M_01060', 'AKUH3B202X_01050',
#               'AKUH3B203J_01080', 'AKUH3B204J_01060', 'AKUH3B205J_01060',
#               'AKUH3B207X_01060', 'AKUH3B208X_01070', 'AKUH4B202J_00980',
#               'AKUH4B204J_01080', 'AKUH4B205J_01060', 'AKUH4B211M_01100',
#               'AKUH4B412M_01140', 'AKUH4B501J_01140', 'AKUH4B502X_01070',
#               'AKUH4B507J_01070', 'AKUH4B507X_01060', 'AKUH4B508X_01060',
#               'MUB42_02675']

# tetR_tags = ['K2W83_RS00575', 'AKUFHON2_01060', #'AAPFHON13_00970', 
#               'AKUG0101_01070', 'AKUG0102_01060', 'AKUG0103_01060', 
#               'AKUG0401_01060', 'AKUG0402_01060', 'AKUG0403_PLPX00280',
#               'AKUG0404_01060', 'AKUG0405_01060', 'AKUG0406_PLPX00290',
#               'AKUG0407_01060', 'AKUG0408_01060', 'AKUG0410_01100',
#               'AKUG0412_01100', 'AKUG0414_01060', 'AKUG0415_01060',
#               'AKUG0417_01090', 'AKUG0420_PLPX00320', 'AKUG0601_01060',
#               'AKUG0602_01060', 'AKUG0702_01060', 'AKUG0801_01060',
#               'AKUG0802_01060', 'AKUG0803_01060', 'AKUG0804_01060',
#               'AKUH1B104J_01060', 'AKUH1B105A_00970', 'AKUH3B101A_01050',
#               'AKUH3B101J_01030', 'AKUH3B102A_01050', 'AKUH3B103J_01050',
#               'AKUH3B103M_PLPX00280', 'AKUH3B104J_01030', 'AKUH3B104X_PLPX00280',
#               'AKUH3B107A_01050', 'AKUH3B109M_01050', 'AKUH3B110M_01050',
#               'AKUH3B111A_PLPX00270', 'AKUH3B111M_01050', 'AKUH3B202X_01040',
#               'AKUH3B203J_01070', 'AKUH3B204J_01050', 'AKUH3B205J_01050',
#               'AKUH3B207X_01050', 'AKUH3B208X_01060', 'AKUH4B202J_00970',
#               'AKUH4B204J_01070', 'AKUH4B205J_01050', 'AKUH4B211M_01090',
#               'AKUH4B412M_01130', 'AKUH4B501J_01130', 'AKUH4B502X_01060',
#               'AKUH4B507J_01060', 'AKUH4B507X_01050', 'AKUH4B508X_01050',
#               'MUB42_02670']
    
    if adhesin == adhesins[0]:
        adh_tags = ['K2W83_RS00570', 'AKUFHON2_01050', #'AAPFHON13_00970', 
                      'AKUG0102_01050', 'AKUG0403_PLPX00290', 'AKUH1B104J_01050', 
                      'AKUH1B105A_00960', 'AKUH3B103M_PLPX00290', 'AKUH3B104J_01020', 
                      'AKUH3B104X_PLPX00300', 'AKUH3B111M_01040', 'AKUH3B202X_01030',
                      'AKUH3B203J_01060', 'AKUH4B202J_00960', 'AKUH4B204J_01060',
                      'AKUH4B412M_01120', 'AKUH4B501J_01120', 'MUB42_02660'] #1 loctag less for HNS-8
        
        pre_adh_tags = ['K2W83_RS00565', 'AKUFHON2_01040', #'AAPFHON13_00970', 
                      'AKUG0102_01040', 'AKUG0403_PLPX00300', 'AKUH1B104J_01040', 
                      'AKUH1B105A_00950', 'AKUH3B103M_PLPX00300', 'AKUH3B104J_01010', 
                      'AKUH3B104X_PLPX00310', 'AKUH3B111M_01030', 'AKUH3B202X_01020',
                      'AKUH3B203J_01050', 'AKUH4B202J_00950', 'AKUH4B204J_01050', 
                      'AKUH4B412M_01110', 'AKUH4B501J_01110', 'MUB42_02655'] #1 loctag less for HNS-8
        # adh_tags = ['K2W83_RS00570', 'AKUFHON2_01050', #'AAPFHON13_00970', 
        #               'AKUG0101_01060', 'AKUG0102_01050', 'AKUG0103_01050', 
        #               'AKUG0401_01050', 'AKUG0402_01050', 'AKUG0403_PLPX00290',
        #               'AKUG0404_01050', 'AKUG0405_01050', 'AKUG0406_PLPX00300',
        #               'AKUG0407_01050', 'AKUG0408_01050', 'AKUG0410_01090',
        #               'AKUG0412_01090', 'AKUG0414_01050', 'AKUG0415_01050',
        #               'AKUG0417_01080', 'AKUG0420_PLPX00330', 'AKUG0601_01050',
        #               'AKUG0602_01050', 'AKUG0702_01050', 'AKUG0801_01050',
        #               'AKUG0802_01050', 'AKUG0803_01050', 'AKUG0804_01050',
        #               'AKUH1B104J_01050', 'AKUH1B105A_00960', 'AKUH3B101A_01040',
        #               'AKUH3B101J_01020', 'AKUH3B102A_01040', 'AKUH3B103J_01040',
        #               'AKUH3B103M_PLPX00290', 'AKUH3B104J_01020', 'AKUH3B104X_PLPX00300',
        #               'AKUH3B107A_01040', 'AKUH3B109M_01040', 'AKUH3B110M_01040',
        #               'AKUH3B111A_PLPX00280', 'AKUH3B111M_01040', 'AKUH3B202X_01030',
        #               'AKUH3B203J_01060', 'AKUH3B204J_01040', 'AKUH3B205J_01040',
        #               'AKUH3B207X_01040', 'AKUH3B208X_01050', 'AKUH4B202J_00960',
        #               'AKUH4B204J_01060', 'AKUH4B205J_01040', 'AKUH4B211M_01080',
        #               'AKUH4B412M_01120', 'AKUH4B501J_01120', 'AKUH4B502X_01050',
        #               'AKUH4B507J_01050', 'AKUH4B507X_01040', 'AKUH4B508X_01040',
        #               'MUB42_02660'] #1 loctag less for HNS-8
        
        # pre_adh_tags = ['K2W83_RS00565', 'AKUFHON2_01040', #'AAPFHON13_00970', 
        #               'AKUG0101_01050', 'AKUG0102_01040', 'AKUG0103_01040', 
        #               'AKUG0401_01040', 'AKUG0402_01040', 'AKUG0403_PLPX00300',
        #               'AKUG0404_01040', 'AKUG0405_01040', 'AKUG0406_PLPX00310',
        #               'AKUG0407_01040', 'AKUG0408_01040', 'AKUG0410_01080',
        #               'AKUG0412_01080', 'AKUG0414_01040', 'AKUG0415_01040',
        #               'AKUG0417_01070', 'AKUG0420_PLPX00340', 'AKUG0601_01040',
        #               'AKUG0602_01040', 'AKUG0702_01040', 'AKUG0801_01040',
        #               'AKUG0802_01040', 'AKUG0803_01040', 'AKUG0804_01040',
        #               'AKUH1B104J_01040', 'AKUH1B105A_00950', 'AKUH3B101A_01030',
        #               'AKUH3B101J_01010', 'AKUH3B102A_01030', 'AKUH3B103J_01030',
        #               'AKUH3B103M_PLPX00300', 'AKUH3B104J_01010', 'AKUH3B104X_PLPX00310',
        #               'AKUH3B107A_01030', 'AKUH3B109M_01030', 'AKUH3B110M_01030',
        #               'AKUH3B111A_PLPX00290', 'AKUH3B111M_01030', 'AKUH3B202X_01020',
        #               'AKUH3B203J_01050', 'AKUH3B204J_01030', 'AKUH3B205J_01030',
        #               'AKUH3B207X_01030', 'AKUH3B208X_01040', 'AKUH4B202J_00950',
        #               'AKUH4B204J_01050', 'AKUH4B205J_01030', 'AKUH4B211M_01070',
        #               'AKUH4B412M_01110', 'AKUH4B501J_01110', 'AKUH4B502X_01040',
        #               'AKUH4B507J_01040', 'AKUH4B507X_01030', 'AKUH4B508X_01030',
        #               'MUB42_02655'] #1 loctag less for HNS-8
        
    elif adhesin == adhesins[1]:
        adh_tags = ['AKUA1003_13820', 'AKUA1202_14880', 'AKUA1401_14120',
                    'AKUA1805_14110', 'K2W83_RS06805', 'AKUFHON2_14900',
                    'AKUG0102_14080', 'AKUG0403_14370', 'VQ058_RS06965',
                    'AKUH1B104J_14360', 'AKUH1B105A_13510', 'AKUH1B302M_14240',
                    'AKUH3B103M_14630', 'AKUH3B104J_14310', 'AKUH3B104X_14630',
                    'AKUH3B111M_13900', 'AKUH3B202X_14170', 'AKUH3B203J_14720',
                    'AKUH3B203M_13780', 'AKUH3B209X_14730', 'AKUH4B202J_14100',
                    'AKUH4B204J_14700', 'AKUH4B402J_13930', 'AKUH4B412M_14600',
                    'AKUH4B501J_14220', 'AKUH4B503X_14010', 'AKUH4B504J_14660',
                    'AKUH4B505J_14210', 'APS55_RS03170'] #Note: only one predicted for G0403 and H3B2-03M? Why?
        
        pre_adh_tags = ['AKUA1003_13830', 'AKUA1202_14890', 'AKUA1401_14130',
                    'AKUA1805_14120', 'K2W83_RS06810', 'AKUFHON2_14910',
                    'AKUG0102_14090', 'AKUG0403_14380', 'VQ058_RS06970',
                    'AKUH1B104J_14370', 'AKUH1B105A_13520', 'AKUH1B302M_14250',
                    'AKUH3B103M_14640', 'AKUH3B104J_14320', 'AKUH3B104X_14640',
                    'AKUH3B111M_13910', 'AKUH3B202X_14180', 'AKUH3B203J_14730',
                    'AKUH3B203M_13790', 'AKUH3B209X_14740', 'AKUH4B202J_14110',
                    'AKUH4B204J_14710', 'AKUH4B402J_13940', 'AKUH4B412M_14610',
                    'AKUH4B501J_14230', 'AKUH4B503X_14020', 'AKUH4B504J_14670',
                    'AKUH4B505J_14220', 'APS55_RS03165']
        
    elif adhesin == adhesins[2]:
        adh_tags = ['AKUA0901_00560', 'AKUA1003_00540', 'AKUA1202_00570',
                    'AKUA1401_00560', 'AKUA1805_00590', 'K2W83_RS00320',
                    'AKUFHON2_00530', 'AKUG0102_00540', 'AKUG0403_00530',
                    'VQ058_RS00290', 'AKUH1B104J_00520', 'AKUH1B105A_00480',
                    'AKUH1B302M_00560', 'AKUH3B103M_00540', 'AKUH3B104J_00520',
                    'AKUH3B104X_00540', 'AKUH3B111M_00510', 'AKUH3B202X_00520',
                    'AKUH3B203J_00530', 'AKUH3B209X_00540', 'AKUH4B111J_00580',
                    'AKUH4B202J_00480', 'AKUH4B204J_00540', 'AKUH4B206J_00580',
                    'AKUH4B402J_00530', 'AKUH4B412M_00590', 'AKUH4B501J_00570',
                    'AKUH4B503X_00500', 'AKUH4B504J_00560', 'AKUH4B505J_00560',
                    'MUB42_02470', 'LDX55_00290', 'APS55_RS02515']
        
        pre_adh_tags = ['AKUA0901_00570', 'AKUA1003_00550', 'AKUA1202_00580',
                        'AKUA1401_00570', 'AKUA1805_00600', 'K2W83_RS00325',
                        'AKUFHON2_00540', 'AKUG0102_00550', 'AKUG0403_00540',
                        'VQ058_RS00295', 'AKUH1B104J_00530', 'AKUH1B105A_00490',
                        'AKUH1B302M_00570', 'AKUH3B103M_00550', 'AKUH3B104J_00530',
                        'AKUH3B104X_00550', 'AKUH3B111M_00520', 'AKUH3B202X_00530',
                        'AKUH3B203J_00540', 'AKUH3B209X_00550', 'AKUH4B111J_00590',
                        'AKUH4B202J_00490', 'AKUH4B204J_00550', 'AKUH4B206J_00590',
                        'AKUH4B402J_00540', 'AKUH4B412M_00600', 'AKUH4B501J_00580',
                        'AKUH4B503X_00510', 'AKUH4B504J_00570', 'AKUH4B505J_00570',
                        'MUB42_02475', 'LDX55_00295', 'APS55_RS02510']

# =============================================================================
# 3. Define paths to inputs and outputs
# =============================================================================
    
    strand_dict = {} #Dictionary to store the strand of the previous CDS
    outfile = f'{outdir}/{adhesin}_intergenic.fna' #Output fna file with the intergenic region
    mafft_outfile = outfile.replace('.fna', '.mafft.fna') #Output alignment
    trimmed_outfile = mafft_outfile.replace('.mafft', '.trimmed.mafft') #Trimmed output alignment
    plot_file = f'{plotdir}/{adhesin}_intergenic_aln.png' #Output alignment plot
    
    get_strains = get_unique_strains(adh_tags) #Retrieve strain names
    
    #Retrieve list of input GenBank and genomic FASTA files
    #1. First, retrieve the paths to the GenBank files (in two different directories)
    infiles = [f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('.gbff') and file.split('_')[0] in get_strains]
    infiles += [f'{indir2}/{file}' for file in os.listdir(indir2) if file.endswith('.gbff') and file.split('_')[0] in get_strains]
    infiles = sorted(infiles) #Sort alphabetically
    #2. Then, retrieve the paths to the FASTA files
    fna_infiles = [file.replace('gbff', 'fna') for file in infiles]
    
# =============================================================================
# 4. Retrieve the positions to save to file
# =============================================================================
    
    #Previous locus tag variable needed to properly retrieve the positions of genes in the reverse strand followed by genes in the forward strand
    prev_loctag = ''
    
    print(f'Retrieving positions for {adhesin}!')
    pos_dict = {} #Dictonary to store the positions
    for file in infiles: #Loop through GenBank files
        with open(file) as handle: #Open input file
            for record in SeqIO.parse(handle, 'genbank'): #Loop through records (contigs) in file
                for cds in record.features: #Loop through features in each record (genes)
                    if cds.type == 'CDS' and 'locus_tag' in list(cds.qualifiers.keys()): #If the gene has an assigned locus tag
                        loctag = cds.qualifiers['locus_tag'][0] #Retrieve the locus tag
                        #If the locus tag corresponds to the gene before the adhesin and the adhesin is in the forward strand
                        if loctag in pre_adh_tags and (cds.location.strand == 1 and 'PLPX' not in loctag):
                            start = int(cds.location.end) #Set the end position of the gene as the start
                            #The strand is considered forward only if the gene is not preceeded by an adhesin in the reverse strand, to set the right positions
                            strand = '+' if not (prev_loctag in adh_tags and strand_dict[prev_loctag] == '-') else '-'
                            start = int(cds.location.end) if strand == '+' else int(cds.location.start) #Set the end position of the gene as the start
                        elif loctag in pre_adh_tags and (cds.location.strand == -1 or 'PLPX' in loctag): #If the gene is in the reverse strand
                            start = int(cds.location.start) #Set the start of the intergenic region to the start of the gene
                            strand = '-'
                        elif loctag in adh_tags and (cds.location.strand == 1 and 'PLPX' not in loctag): #If the locus tag is the adhesin and it is in the forward strand
                            end = int(cds.location.start) + 12 #Set the end of the segment to the gene start location + 12 nucleotides
                            strand = '+'
                        elif loctag in adh_tags and (cds.location.strand == -1 or 'PLPX' in loctag): #If the locus tag corresponds to the adhesin and it is in the reverse strand
                            end = int(cds.location.end) - 12 #Set the end of the segment to the end of the gene - 12 nucleotides
                            strand = '-'
                        prev_loctag = loctag #Store previous locus tag
                        if loctag in adh_tags: #If the locus tag corresponds to an adhesin
                            strand_dict[loctag] = strand #Store strand
        strain = os.path.basename(file).split('_')[0] #Retrieve the strain name
        pos_dict[strain] = (start, end, strand) #Save the position
        print(f'{strain}: {start}-{end} ({strand}), {abs(end-start)} nucleotides.')
        
# =============================================================================
# 5. Save positions of interest to file
# =============================================================================
        
    with open(outfile, 'w') as intergenic: #Open output file
        for fna in fna_infiles: #Loop through fna files
            strain = os.path.basename(fna).split('_')[0] #Retrieve strain name
            if adhesin == adhesins[0] and pos_dict[strain][1] < pos_dict[strain][0]: #If the start position is bigger than the end position (reverse strand)
                plasmid = True #Set boolean to false
            else: plasmid = False #Else, set boolean to false
            with open(fna) as handle: #Open FASTA file
                for record in SeqIO.parse(handle, 'fasta'): #Loop through records in the file
                    if not plasmid: #If the gene is not in the plasmid
                        start = min(pos_dict[strain][0], pos_dict[strain][1]) #Get lowest value
                        end = max(pos_dict[strain][0], pos_dict[strain][1]) #Get highest value
                        fna_seq = record.seq[start:end] #Retrieve segment
                        if pos_dict[strain][2] == '-': #If the gene is in the reverse strand
                            fna_seq = fna_seq.reverse_complement() #Get reverse complement of segments in the reverse strand
                        new_record = SeqRecord(fna_seq, id = f'{strain}', description = '') #Create new record with the sequence
                        intergenic.write(as_fasta(new_record)) #Write as fasta
                        break

                
# =============================================================================
# 6. Align the sequences
# =============================================================================
                
    #Align the intergenic regions and trim the alignments to remove gappy positions
    subprocess.run(f'mafft-linsi {outfile} > {mafft_outfile}', shell = True)
    subprocess.run(f'trimal -in {mafft_outfile} -out {trimmed_outfile} -gappyout -fasta', shell = True)

# =============================================================================
# 7. Plot the alignments and write the consensus to file
# =============================================================================
    
    #Create the object to be plotted and assign colors to the plot
    mv = MsaViz(trimmed_outfile, wrap_length = 100, color_scheme = 'Identity', 
                show_consensus = True, consensus_color = '#A4BF19')
    mv.savefig(plot_file, dpi = 300) #Save the plot to a file
    
    consensus = mv._consensus_seq.upper() #Retrieve the consensus sequence
    #Convert it to a Biopython record
    consensus_record = SeqRecord(Seq(consensus), id = adhesin,
                                 description = 'intergenic consensus')

    with open(consensus_file, 'a') as handle: #Open output file with consensus sequences
        handle.write(as_fasta(consensus_record)) #Write the consensus to the file
        
# =============================================================================
# 8. Align the consensus sequences and plot the alignment
# =============================================================================

#Run MAFFT to align all the consensus to each other
subprocess.run(f'mafft-linsi {consensus_file} > {mafft_consensus}', shell = True)

#Plot the consensus
mv = MsaViz(mafft_consensus, wrap_length = 100, color_scheme = 'Identity', 
            show_consensus = True, consensus_color = '#A4BF19')
mv.savefig(consensus_plot, dpi = 300)