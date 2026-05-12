#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nay 12, based on script from Mon Sep  4 10:28:42 2023

This is the code to generate a Blast comparison plot in PygenomeViz focusing
on the regions where the most differentially expressed adhesins are located.

@author: Marina Mota-Merlo
"""

# =============================================================================
# Import required modules
# =============================================================================

from matplotlib.lines import Line2D
from Bio.SeqFeature import SimpleLocation
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank as gbk_read
from Bio import GenBank as gbk
import os
import subprocess
import pandas as pd

# =============================================================================
# Color dictionary for strain names
# =============================================================================
leaf_color = {'A0901': '#D55E00', 'A1001': '#771853', 'A1003': '#0072B2', 
              'A1202': '#33B18F', 'A1401': '#33B18F', 'A1404': '#FF74D6', 
              'A1805': '#33B18F', 'G0101': '#0072B2', 'G0403': '#33B18F', 
              'Fhon2': '#0072B2', 'H1B104J': '#0072B2', 'H1B105A': '#0072B2', 
              'H1B302M': '#0072B2', 'H3B104J': '#0072B2', 'H3B104X': '#0072B2',
              'H3B101A': '#0072B2', 'H3B202X': '#0072B2', 'H3B203J': '#0072B2', 
              'H3B203M': '#D55E00', 'H3B206M': '#D55E00', 'H3B209X': '#33B18F', 
              'H4B111J': '#D55E00', 'H4B202J': '#0072B2', 'H4B204J': '#0072B2', 
              'H4B206J': '#D55E00', 'H4B402J': '#0072B2', 'H4B405J': '#D55E00',
              'H4B406M': '#D55E00', 'H4B412M': '#0072B2', 'H4B501J': '#0072B2', 
              'H4B503X': '#0072B2', 'H4B504J': '#33B18F', 'H4B505J': '#33B18F', 
              'MP2': '#33B18F', 'IBH001': '#D55E00', 'DSMZ12361': '#0072B2', 
              'HNS8': 'black', 'GYUN333': 'black', 'Fhon13': '#79443B'}

# =============================================================================
# Locus tags
# =============================================================================
#Note: All representative strains will be included in the final plot

adh1 = ['K2W83_RS00570', 'AKUFHON2_01050', 'AKUG0102_01050', 
        'AKUH1B104J_01050', 'AKUH1B105A_00960', 'AKUH3B101A_01040', 
        'AKUH3B104J_01020', 'AKUH3B202X_01030', 'AKUH3B203J_01060', 
        'AKUH4B202J_00960', 'AKUH4B204J_01060', 'AKUH4B412M_01120', 
        'AKUH4B501J_01120', 'MUB42_02665']
#Representative strains
#Use the "product" section of the GeneBank, check that it contains NrdI

adh2 = ['AAPFHON13_01040', 'AKUA1805_01330', 'K2W83_RS00655',
        'AKUH3B101A_01250', 'AKUH3B202X_01200', 'AKUH3B202X_01250', 
        'AKUH3B203J_01270', 'AKUH3B203M_01290', 'AKUH4B504J_01250'] #Representative strains with the adhesin

collagen = ['AKUA0901_00550', 'AKUA1003_00530', 'AKUA1202_00540',
        'AKUA1401_00550', 'AKUA1805_00580', 'K2W83_RS00315',
        'AKUFHON2_00520', 'AKUG0101_00540', 'AKUG0403_00520',
        'VQ058_RS00285', 'AKUH1B104J_00510', 'AKUH1B105A_00470',
        'AKUH1B302M_00550', 'AKUH3B103M_00530', 'AKUH3B104J_00510',
        'AKUH3B104X_00530', 'AKUH3B101A_00500', 'AKUH3B202X_00510',
        'AKUH3B203J_00520', 'AKUH3B209X_00530', 'AKUH4B111J_00570',
        'AKUH4B202J_00470', 'AKUH4B204J_00530', 'AKUH4B206J_00570',
        'AKUH4B402J_00520', 'AKUH4B412M_00580', 'AKUH4B501J_00560',
        'AKUH4B503X_00490', 'AKUH4B504J_00550', 'AKUH4B505J_00550',
        'MUB42_02465', 'LDX55_00285', 'APS55_RS02520']
#Strains that have both adhesins + representatives
#Use the "product" section of the GeneBank, check that it contains NrdI

gtf2 = ['AKUA0901_00560', 'AKUA1003_00540', 'AKUA1202_00570',
        'AKUA1401_00560', 'AKUA1805_00590', 'K2W83_RS00320',
        'AKUFHON2_00530', 'AKUG0101_00550', 'AKUG0403_00530',
        'VQ058_RS00290', 'AKUH1B104J_00520', 'AKUH1B105A_00480',
        'AKUH1B302M_00560', 'AKUH3B103M_00540', 'AKUH3B104J_00520',
        'AKUH3B104X_00540', 'AKUH3B101A_00510', 'AKUH3B202X_00520',
        'AKUH3B203J_00530', 'AKUH3B209X_00540', 'AKUH4B111J_00580',
        'AKUH4B202J_00480', 'AKUH4B204J_00540', 'AKUH4B206J_00580',
        'AKUH4B402J_00530', 'AKUH4B412M_00590', 'AKUH4B501J_00570',
        'AKUH4B503X_00500', 'AKUH4B504J_00560', 'AKUH4B505J_00560',
        'MUB42_02470', 'LDX55_00290', 'APS55_RS02515']

adh3 = ['AKUA1003_13810', 'AKUA1202_14870', 'AKUA1401_14110',
            'AKUA1805_14100', 'K2W83_RS06800', 'AKUFHON2_14890',
            'AKUG0101_14160', 'AKUG0403_14360', 'VQ058_RS06960',
            'AKUH1B104J_14350', 'AKUH1B105A_13500', 'AKUH1B302M_14230',
            'AKUH3B103M_14620', 'AKUH3B104J_14300', 'AKUH3B104X_14620',
            'AKUH3B101A_14570', 'AKUH3B202X_14160', 'AKUH3B203J_14710',
            'AKUH3B203M_13770', 'AKUH3B209X_14720', 'AKUH4B202J_14090',
            'AKUH4B204J_14690', 'AKUH4B402J_13920', 'AKUH4B405J_14560', 
            'AKUH4B406M_14860', 'AKUH4B412M_14590', 'AKUH4B501J_14210', 
            'AKUH4B503X_14000', 'AKUH4B504J_14650', 'AKUH4B505J_14200', 
            'APS55_RS03170'] #Note: only one predicted for G0403 and H3B2-03M? Why?

adh4 = ['AKUA1003_13820', 'AKUA1202_14880', 'AKUA1401_14120',
            'AKUA1805_14110', 'K2W83_RS06805', 'AKUFHON2_14900',
            'AKUG0101_14170', 'AKUG0403_14370', 'VQ058_RS06965',
            'AKUH1B104J_14360', 'AKUH1B105A_13510', 'AKUH1B302M_14240',
            'AKUH3B103M_14630', 'AKUH3B104J_14310', 'AKUH3B104X_14630',
            'AKUH3B101A_14580', 'AKUH3B202X_14170', 'AKUH3B203J_14720',
            'AKUH3B203M_13780', 'AKUH3B209X_14730', 'AKUH4B202J_14100',
            'AKUH4B204J_14700', 'AKUH4B402J_13930', 'AKUH4B412M_14600',
            'AKUH4B501J_14220', 'AKUH4B503X_14010', 'AKUH4B504J_14660',
            'AKUH4B505J_14210', 'APS55_RS03175'] 

adh5 = ['AKUA0901_14630', 'AKUA1001_13630', 'AKUH3B206M_14150', 'AKUH4B111J_14900',
        'AKUH4B206J_14800', 'LDX55_06980']

#Gene order when plotting
phylo_order = {1: 'H4B4-12M', 2: 'H4B5-01J', 3: 'G0101', 4: 'H1B1-04J',
               5: 'H4B5-03X', 6: 'H4B4-02J', 7: 'H1B3-02M', 8: 'H3B2-02X',
               9: 'H3B1-04X', 10: 'A1003', 11: 'H3B1-01A', 12: 'H3B2-03J',
               13: 'H4B2-04J', 14: 'Fhon2', 15: 'H4B2-02J', 16: 'H1B1-05A',
               17: 'DSMZ12361', 18: 'H3B1-04J', 19: 'H3B2-06M', 20: 'H4B4-05J',
               21: 'H4B1-11J', 22: 'H4B4-06M', 23: 'H3B2-03M', 24: 'A0901',
               25: 'H4B2-06J', 26: 'IBH001', 27: 'H3B2-09X', 28: 'H4B5-05J',
               29: 'A1401', 30: 'A1202', 31: 'H4B5-04J', 32: 'MP2', 33: 'A1805',
               34: 'G0403', 35: 'A1404', 36: 'A1001', 37: 'GYUN-333', 38: 'HNS-8'}

#Colors for the genes of interest
color_dict = {'adh1': '#BBE36A', 'adh2': '#5FB477', 'gtf2': '#72A3E0', 
              'collagen': '#E072A3', 'TF': '#5603AD', 'FK': '#E0777D',
              'SK': '#B973FF', 'RE': '#FF802B', 'large': '#037971', 
              'adh3': '#F7ED5C', 'adh4': '#F5BB64', 'adh5': '#F59264'}

# =============================================================================
# In this section, we run Blast between pairs of strains following the 
# phylogeny.
# =============================================================================

adhesin_regions = ['MucBP+LPXTG', 'collagen_binding_Gtf2', 'MubB2+LPXTG'] #List of adhesin regions to retrieve
ref_genes = ['NrdI', 'WalK', 'omoserine kinase'] #List of reference genes that are present across all strains near the region
projdir = os.path.expanduser('~') + '/mucoid_project/adhesins' #Working directory
outpath = f'{projdir}/blast_tabs' #Output path to save BLAST comparisons
folder_path = os.path.expanduser('~') + '/Akunkeei_files/fna' #Path to the input FNA files
folder_path2 = os.path.expanduser('~') + '/Akunkeei_files/gbff' #Path to the input GBFF files

fna_files = [f'{folder_path}/{strain}_genomic.fna' for strain in phylo_order.values()] #Retrieve the paths to all FNA files

#Modify the paths of a couple of strains for which the files are in a different directory
fna_files[-2] = fna_files[-2].replace('/fna', '/new_genomes/fna')
fna_files[-1] = fna_files[-1].replace('/fna', '/new_genomes/fna')

length_dict = {} #Dictionary to store chromosome lengths

def get_tabs(fnas, outpath):
    """
    Function to run BLASTn between pairs of strains
    """
    if not os.path.exists(outpath): #If the output directory does not exist
        os.makedirs(outpath) #Create it
    for i in range(1, len(fnas)): #Loop through all the input files
        strain1_file = fnas[i-1] #Retrieve the first file in the comparison
        strain2_file = fnas[i] #Retrieve the second file
        outfile = f'{outpath}/{i}.tab' #Output file
        if i < 10: #If the comparison number is below 10
            outfile = outfile.replace(f'{i}', f'0{i}') #Add a 0 before (for sorting purposes)
        
        #Run BLAST to calculate matches between strains
        subprocess.run(f'blastn -query {strain2_file} -subject {strain1_file} -outfmt 7 -out {outfile}', shell = True)
        
get_tabs(fna_files, outpath)

def get_info(fnas, folder_path2, gene_name, adhesin_region, length = 0):
    '''
    Function to get the start and end coordinates to be plotted.
    
    Input parameters:
        @param phylo_dict (dict, {str: int}): Dictionary indicating the order 
        in which the strains should be plotted. This parameter should 
        correspond to the order for the Blast comparisons.
        @param folder_path2 (str): Path where the input GenBank files are.
        @param length (int): Length of the segment to be plotted.
    
    Output parameters:
        @param accession_strain (dict, {str: [str, ...]}): Dictionary with 
        strain names as keys and NCBI accessions as values.
        @param pos_dict (dict, {str: (int, int)}): Dictionary with tuples
        indicating start and end positions.
        @lengths_dict (dict, {str: int}): Dictionary with the names of the strains
        as keys and the lengths of the chromosome as values.
    '''

    accession_strain = {} #Temporary variable to store NCBI accessions
    pos_dict = {} #Dictionary to store the total genome length of each strain
    lengths_dict = {} #Diccionary of start and end positions
    for fna in fnas: #Loop through strain names
        genbank_file = fna.replace('fna', 'gbff') #Set name of GenBank file
        strain = os.path.basename(genbank_file).replace('_genomic.gbff', '')
        with open(genbank_file) as handle: #Read GenBank file
            for record in gbk.parse(handle): #Loop through the records (chromosomes or plasmids) in the file
                if strain not in accession_strain.keys(): #If the strain is not in the accession dictionary (we only get the first record)
                    accession_strain[strain] = f'{record.accession[0]}.1' #Associate accession to strain name
                    genome_length = len(record.sequence) #Get length of the record
                    lengths_dict[strain] = genome_length #Initialize length dictionary
                    
                # if strain not in ['IBH001', 'MP2', 'DSMZ12361']: #Retrieve start position in the strains
                for feature in record.features:  #Loop through features (mostly CDS) in the record
                    # gene_name = [qual for qual in feature.qualifiers if 'gene' in qual.key] #Get the four-letter gene name
                    product = [qual for qual in feature.qualifiers if 'product' in qual.key]
                    if len(product) > 0 and gene_name in product[0].value:
                        start = int(feature.location.split('..')[0].replace('complement(', '')) #Set start location to the beginning of the reference gene
                    # if len(gene_name) > 0: #If there is at least one qualifier
                    #     gene = gene_name[0] #Set the gene name to the first qualifier (gene name or locus tag)
                        if adhesin_region == adhesin_regions[2] and strain == 'HNS-8':
                            break

            if (strain != 'MP2' and adhesin_region != adhesin_regions[0]) or (strain == 'MP2' and adhesin_region == adhesin_regions[0]): #If the strain is not MP2 (assembly of the reverse strand)
                start = lengths_dict[strain] - (start + length) #Set the start to the opposite strand
                end = start + length #Set the end to the start + segment length
                pos_dict[strain] = (start, end) #Save the new position to a dictionary
                if strain == 'MP2':
                    start = end + length
                    pos_dict[strain] = (end, start)
            else: #If the strain is MP2
                end = start + length #Calculate the beginning point of the segment
                pos_dict[strain] = (start, end) #Store the results in a dictionary
                if strain == 'MP2':
                    end = start-length
                    pos_dict[strain] = (end, start)
    return accession_strain, pos_dict, lengths_dict

for i in range(len(adhesin_regions)):
    
    true_dict = {gene: False for gene in color_dict.keys()}
    
    adhesin_region = adhesin_regions[i]
    ref = ref_genes[i]
    outfig = f'{projdir}/plots/pyGenomeViz/{adhesin_region}_region.svg'
    
    if not os.path.exists(os.path.dirname(outfig)):
        os.makedirs(os.path.dirname(outfig))
    
    # =============================================================================
    # In this section, we establish the start and end positions for every strain
    # and save them to a dictionary. I also get the chromosome IDs to change them
    # to strain IDs.
    # =============================================================================
        
    acc, pos, lengths = get_info(fna_files, folder_path2, ref, adhesin_region, 70000)
    
    
    # =============================================================================
    # Now, ready for the plotting!
    # =============================================================================
    
    # Set plot style
    gv = GenomeViz(
        fig_track_height = 0.42, #Height of the tracks with the representation of the CDS
        link_track_ratio = 0.8, #Size ratio between the links (Blast comparisons) and the tracks
        track_align_type = 'center', #Align tracks to the center
        )
    
    gv.set_scale_bar(scale_size_label=(5000, '5 Kb')) #Set a legend for the scale of the plot (5kb bar)
    
    for fna in fna_files: #Loop through strain names in the order dictionary
        strain = os.path.basename(fna).replace('_genomic.fna', '') #Get the name of a strain
        genbank_file = fna.replace('fna', 'gbff') #Get GenBank file based on strain name
        genbk = gbk_read(genbank_file) #Read GenBank file with PyGenomeViz
        segments = dict(region1=(pos[strain][0], pos[strain][1])) #Retrieve the segment to be plotted
    
        track = gv.add_feature_track(name = genbk.name.replace('_genomic', '').replace('12361', ''), #Create track (use DSMZ as strain name for DSMZ12361)
                                     segments = segments, #Add segments to track
                                     label_kws = dict(color = leaf_color[strain.replace('-', '')])) #Add strain name color based on phylogroup
        chromosome_length = len(genbk.records[0].seq)
        length_dict[strain] = chromosome_length
        for segment in track.segments: #Loop through segments in the track
            if (strain != 'MP2' and adhesin_region != adhesin_regions[0]) or (strain == 'MP2' and adhesin_region == adhesin_regions[0]): #If the strain is not MP2
                target_range = (lengths[strain] - segment.range[1], 
                                lengths[strain] - segment.range[0]) #Get the target region on the opposite strand
                features = genbk.extract_features(feature_type = 'CDS', 
                                                  target_range = target_range) #Extract features from GenBank file
    
                for feature in features: #Loop through features
                    new_end = lengths[strain] - feature.location.start #Reverse start
                    new_start = lengths[strain] - feature.location.end #Reverse end
                    print(new_start, new_end)
                    new_strand = feature.location.strand*-1 #Reverse strand
                    feature.location = SimpleLocation(new_start, new_end, 
                                                      new_strand) #Update feature position
            else: #If the strain is MP2
                features = genbk.extract_features(feature_type = 'CDS', 
                                              target_range = segment.range)  #Extract the features directly, based on the segment range
            #Loop through CDS
            for cds in features: #Loop through CDS
                protstart = int(cds.location.start) #Get CDS start
                end = int(cds.location.end) #Get CDS end
                strand = cds.location.strand #Get strand
                color = '#E3DAC9' #Set color of most CDS
                gene_name  = '' #Initialize gene name
                if end - protstart > 7000 or cds.qualifiers['locus_tag'][0] in ['AKUH1B104J_01280', 'AKUH1B104J_01290', 'AKUA1805_01370', 'AKUA1805_01380', 'AKUA2101_01370', 'AKUA2101_01380', 'AKUH1B105A_01180']:
                    color = color_dict['large']
                    true_dict['large'] = True
                # elif 'glycosyl hydrolase' in cds.qualifiers['product'][0] and cds.qualifiers['locus_tag'][0] not in gtf2: #Retrieve the set of genes that were manually annotated in the GenBanks
                #     color = color_dict['GH']
                elif 'fructokinase' in cds.qualifiers['product'][0].lower() or '6-phosphate' in cds.qualifiers['product'][0] or 'ROK' in cds.qualifiers['product'][0]:
                    color = color_dict['FK']
                    true_dict['FK'] = True
                elif 'kinase' in cds.qualifiers['product'][0].lower():
                    color = color_dict['SK']
                    true_dict['SK'] = True
                elif 'transcription factor' in cds.qualifiers['product'][0] or 'regulator' in cds.qualifiers['product'][0] or 'sensory' in cds.qualifiers['product'][0].lower():
                    color = color_dict['TF']
                    true_dict['TF'] = True
                elif cds.qualifiers['locus_tag'][0] in collagen:
                    color = color_dict['collagen']
                    true_dict['collagen'] = True
                elif cds.qualifiers['locus_tag'][0] in gtf2:
                    color = color_dict['gtf2']
                elif cds.qualifiers['locus_tag'][0] in adh1:
                    color = color_dict['adh1']
                    true_dict['adh1'] = True
                elif cds.qualifiers['locus_tag'][0] in adh2:
                    color = color_dict['adh2']
                    true_dict['adh2'] = True
                    if cds.qualifiers['locus_tag'][0] == 'AKUH3B202X_01250':
                        color = '#8DCC70'
                elif cds.qualifiers['locus_tag'][0] in adh3:
                    color = color_dict['adh3']
                    true_dict['adh3'] = True
                elif cds.qualifiers['locus_tag'][0] in adh4:
                    color = color_dict['adh4']
                elif cds.qualifiers['locus_tag'][0] in adh5:
                    color = color_dict['adh5']
                elif 'restriction' in cds.qualifiers['product'][0] or 'nuclease' in cds.qualifiers['product'][0]:
                    color = color_dict['RE']
                if strain == 'Fhon13' and 'product' in cds.qualifiers.keys():
                    possible_name = cds.qualifiers['product'][0].split(' ')[0]
                    if len(cds.qualifiers['product'][0]) < 7 or (len(possible_name) < 7 and not (possible_name[0].isupper() and possible_name[1].islower())):
                        gene_name = possible_name
                if 'gene' in cds.qualifiers.keys():
                    gene_name = cds.qualifiers['gene'][0] #When possible, add four-letter gene name
                    if 'transposase' not in cds.qualifiers['product'][0]: #If the gene is not a transposon
                        gene_name = cds.qualifiers['gene'][0]
                    else: #If the gene is a transposon
                        color = 'black' #Color it in black
                        
                elif 'transposase' in cds.qualifiers['product'][0]: #Color transposases in black
                    gene_name = ''
                    color = 'black'
                 
                if segment.start <= protstart <= end <= segment.end: #If the CDS is inside the segment to be plotted
                    segment.add_feature(protstart, end, strand, label = gene_name, #Add CDS to segment (position and label)
                                      plotstyle = 'bigarrow', fc = color, lw = 1, arrow_shaft_ratio = 1, #Set arrow style, CDS color, CDS line width and arrow vs shaft ratio
                                      text_kws = dict(color = 'black', rotation = 45, #Set label properties (label and rotation)
                                                      size = 15, ymargin = 0, #Set label properties (text size and distance from the CDS)
                                                      vpos = 'top', hpos = 'left')) #Set label properties (vertical and horizontal position)
                
            if (strain != 'MP2' and adhesin_region != adhesin_regions[0]) or (strain == 'MP2' and adhesin_region == adhesin_regions[0]): #To change the direction of the plot depending on the strand
                segment.add_sublabel(f'{chromosome_length-segment.start:,} - {chromosome_length-segment.end:,} bp')  #Add text indicating segment range to the plot
            else:
                segment.add_sublabel(f'{segment.start:,} - {segment.end:,} bp')  #Add text indicating segment range to the plot
    
            track.align_label = True #Align track label (strain name) to track
            track.set_segment_sep() #Set separator (//) between segments
    
            
    # # =============================================================================
    # # Here I modify sligthly the tab files that I created.
    # # =============================================================================
    for i in range(1, len(phylo_order.keys())): #Loop through strains in the order in which they will be plotted
        if i < 10:
            tab = f'{outpath}/0{i}.tab'
        else:
            tab = f'{outpath}/{i}.tab'
        with open(tab) as tabfile: #Open Blast comparison file for each strain
            k = 0 #Initialize parameter to retrieve column names
            for line in tabfile: #Loop through lines in tab file
                    k += 1 #Increase the value of k
                    #Generate a list with the fields:
                    if 'Fields' in line: #If the line contains the string Fields
                        headers = line #Assign the header string to a variable
                    elif k > 10: #If we pass the line
                        break #End loop
        
        header_list = headers.split(',') #Divide comma-separated fields in the header string to create a list
        header_list[0] = header_list[0].replace('# Fields: ', '') #Remove the start of the line
        header_list = [header.strip().strip('\n') for header in header_list] #Remove spaces and line breaks from beginning
        
        #Read Blast comparisons as tab-separated file, skip 5 rows and don't set the remaining rows as headers
        file_df = pd.read_csv(tab, sep = '\t', skiprows = 5, header = None)
        file_df = file_df[:-1] #Remove the last row, that also contains a comment
        file_df.columns = header_list #Add column names
        strain2 = file_df.values[0, 0] #Retrieve the accession of the query strain
        strain2_name = list(acc.keys())[list(acc.values()).index(strain2)] #Using the accession, retrieve the key (strain name)
        file_df = file_df.replace(strain2, strain2_name) #Replace the accession with the strain name in the dataframe
        strain1 = file_df.values[0, 1] #Do the same for the subject strain
        strain1_name = list(acc.keys())[list(acc.values()).index(strain1)]
        file_df = file_df.replace(strain1, strain1_name)
        print(f'Loading comparison between {strain1_name} and {strain2_name}') #Print information about the comparison being made
        
        query_pos = pos[phylo_order[i + 1]] #Retrieve start and end position of the query segment
        subject_pos = pos[phylo_order[i]] #Retrieve start and end position of the subject segment
        if adhesin_region != adhesin_regions[0]: #Boundary positions have to be defined differently for the reverse strand
            query_pos = (length_dict[strain2_name] - query_pos[1], length_dict[strain2_name] - query_pos[0]) #Retrieve start and end position of the query segment
            subject_pos = (length_dict[strain1_name] - subject_pos[1], length_dict[strain1_name] - subject_pos[0]) #Retrieve start and end position of the subject segment
        
        tab_df = file_df.copy() #Copy file_df to a new dataframe
        subject_check = False #Check that states whether the subject has one or two segments to plot
        query_check = False #Check that states whether the query has one or two segments to plot
        
        # Figure out how to loop through segments to drop genes (idea: compare strain name with tab df and retrieve those tracks, then loop through segments)
        tab_df = tab_df.drop(tab_df.index[tab_df['q. end'] < query_pos[0]]) #Drop CDS that start after the end of the query segment
        tab_df = tab_df.drop(tab_df.index[tab_df['q. start'] > query_pos[1]]) #Drop CDS that end before the start of the query segment
        tab_df.loc[tab_df['q. start'] < query_pos[0], 'q. start'] = query_pos[0] #Trim the start of matches that overlap with the segment, but start before the segment
        tab_df.loc[tab_df['q. end'] > query_pos[1], 'q. end'] = query_pos[1] #Trim the end of matches that overlap with the segment, but end after the segment
        
        #Do the same for the subject as for the query
        tab_df = tab_df.drop(tab_df.index[tab_df['s. end'] < subject_pos[0]])
        tab_df = tab_df.drop(tab_df.index[tab_df['s. start'] > subject_pos[1]])
        tab_df.loc[tab_df['s. end'] > subject_pos[1], 's. end'] = subject_pos[1]
        tab_df.loc[tab_df['s. start'] < subject_pos[0], 's. start'] = subject_pos[0]
        
        tab_df = tab_df[~tab_df['query acc.ver'].str.startswith('#')] #Remove rows with query strain names starting with #
        tab_df = tab_df.reset_index() #Reset the index of the dataframe
        
        for j in range(len(tab_df)):
            #Get strain names from dataframe and shorten name of DSMZ
            strain1_simp = tab_df.loc[j, 'query acc.ver'].replace('12361', '')
            strain2_simp = tab_df.loc[j, 'subject acc.ver'].replace('12361', '')
            
            #Define the percentage of identity (to color the links)
            identity = tab_df.loc[j, '% identity']
            
            #Define the links (strain name, segment, start, end)
            if adhesin_region != adhesin_regions[0]: #Link positions have to be defined differently if located in the reverse strand
                link1 = (strain1_simp, 'region1', 
                         -1*(tab_df.loc[j, 'q. start'] - length_dict[strain2_name]), 
                         -1*(tab_df.loc[j, 'q. end'] - length_dict[strain2_name]))
                link2 = (strain2_simp, 'region1', 
                         -1*(tab_df.loc[j, 's. start'] - length_dict[strain1_name]),
                         -1*(tab_df.loc[j, 's. end'] - length_dict[strain1_name]))
            else:
                link1 = (strain1_simp, 'region1', tab_df.loc[j, 'q. start'], tab_df.loc[j, 'q. end'])
                link2 = (strain2_simp, 'region1', tab_df.loc[j, 's. start'], tab_df.loc[j, 's. end'])
            
            #Add the links to the gv plot, especifying the colors of Blast matches (v indicates the variable setting color, and vmin is the minimum value)
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = 60, curve = True, alpha = 0.7) #Curve makes the matches form curves
                
    #Set legend for the Blast matches (two adjacent colorbars that show the colors of forward and reverse matches with a minimum of 50% identity)        
    gv.set_colorbar(['grey', 'red'], vmin = 60, bar_height = 0.05, 
                    tick_labelsize = 16, alpha = 0.7)
    fig = gv.plotfig(dpi = 400)
    
    #Create legend text (label) and icons (marker sets the size, color sets marker color, ms sets marker size and ls sets border size)
    #The first two variables are the data to be plotted in the x and y axis, and therefore are empty
    handles = [
        Line2D([], [], marker="", color='black', label="Tracks", ms=20, ls="none"),
        Line2D([], [], marker=">", color='#E3DAC9', label="CDS", ms=20, ls="none")
        ]
    
    #Assign legend depending on which genes are found in the region
    if true_dict['adh1'] == True:
        handles += [Line2D([], [], marker=">", color=color_dict['adh1'], label="Adhesin, LPXTG 8a", ms=20, ls="none")]
    
    if true_dict['adh2'] == True:
        handles += [
            Line2D([], [], marker=">", color='#8DCC70', label="Adhesin, ambiguous", ms=20, ls="none"),
            Line2D([], [], marker=">", color=color_dict['adh2'], label="Adhesin, LPXTG 8b", ms=20, ls="none"),
            Line2D([], [], marker=">", color=color_dict['large'], label="SH3b-containing surface protein", ms=20, ls="none")
            ]
        
    elif true_dict['collagen'] == True:
        handles += [
            Line2D([], [], marker=">", color=color_dict['adh1'], label="Adhesin, LPXTG 8a", ms=20, ls="none"),
            Line2D([], [], marker=">", color=color_dict['collagen'], label="Adhesin, collagen-binding", ms=20, ls="none"),
            Line2D([], [], marker=">", color=color_dict['gtf2'], label="Glycosyl hydrolase, family 2", ms=20, ls="none")
            ]
        
    elif true_dict['adh3'] == True:
        handles += [
            Line2D([], [], marker=">", color=color_dict['adh3'], label="Adhesin, LPXTG 3", ms=20, ls="none"),
            Line2D([], [], marker=">", color=color_dict['adh4'], label="Adhesin, LPXTG 4", ms=20, ls="none"),
            Line2D([], [], marker=">", color=color_dict['adh5'], label="Adhesin, LPXTG 5", ms=20, ls="none")
            ]
        
    if true_dict['TF'] == True:
        handles += [Line2D([], [], marker=">", color=color_dict['TF'], label="Transcription factor", ms=20, ls="none")]
        
    if true_dict['SK'] == True:
        handles += [Line2D([], [], marker=">", color=color_dict['SK'], label="Sensory kinase", ms=20, ls="none")]
        
    if true_dict['FK'] == True:
        handles += [Line2D([], [], marker=">", color=color_dict['FK'], label="Genes for sugar metabolism", ms=20, ls="none")]
    
    handles += [ #Add the basic legend (strain phylogroup colors and common genes)
        Line2D([], [], marker=">", color=color_dict['RE'], label="Restriction endonuclease", ms=20, ls="none"),
        Line2D([], [], marker=">", color='black', label="Transposase", ms=20, ls="none"),
        Line2D([], [], marker="", color='#771853', label="", ms=20, ls="none"),
        Line2D([], [], marker="", color='black', label="Matches", ms=20, ls="none"),
        Line2D([], [], marker="s", color='grey', label="Forward match", ms=20, ls="none"),
        Line2D([], [], marker="s", color='red', label="Reverse match", ms=20, ls="none"),
        Line2D([], [], marker="", color='#771853', label="", ms=20, ls="none"),
        Line2D([], [], marker="", color='#771853', label="Strain names", lw=2, ms=20, ls="none"),
        Line2D([], [], marker="X", color='#0072B2', label="Phylogroup A", ms=15, ls="none"),
        Line2D([], [], marker="X", color='#33B18F', label="Phylogroup B", ms=15, ls="none"),
        Line2D([], [], marker="X", color='#D55E00', label="Phylogroup C", ms=15, ls="none"),
        Line2D([], [], marker="X", color='#FF74D6', label="Phylogroup E", ms=15, ls="none"),
        Line2D([], [], marker="X", color='#771853', label="Phylogroup F", ms=15, ls="none"),
        Line2D([], [], marker="X", color='black', label="Unknown phylogroup", ms=15, ls="none"),
        # Line2D([], [], marker="X", color='#79443B', label="A. apinorum", ms=15, ls="none")
        ]
    
    #Assign legent to figure. bbox_to_anchor sets the position, frameon removes the 
    #frame (border) of the legend box, and labelspacing increases vertical space between legends
    legend = fig.legend(handles=handles, bbox_to_anchor=(1.35, 1), frameon = False,
                        fontsize = 16)
        
    fig.savefig(outfig) #Save figure to SVG
    fig.savefig(outfig.replace('svg', 'png')) #Save figure to PNG
    fig.savefig(outfig.replace('svg', 'pdf')) #Save figure to PDF
    fig.savefig(outfig.replace('svg', 'tiff')) #Save figure to TIFF
    gv.savefig_html(outfig.replace('svg', 'html')) #Save figure to HTML
