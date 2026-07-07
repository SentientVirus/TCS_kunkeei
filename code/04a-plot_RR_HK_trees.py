#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 13:52:47 2025, modified from a script created on Mon Dec 
20 11:06:54 2021 by Marina Mota-Merlo.

This script creates a plot of the adhesin phylogeny, using midpoint rooting.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

from ete4 import Tree
from ete4.treeview import TreeStyle, NodeStyle, TextFace
import os
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
log = f'{logdir}/04a-plot_RR-TF_HK_trees.log'
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

#Dictionary to color strains by phylogroup
leaf_color = {'A0901': '#D55E00', 'A1001': '#771853', 'A1002': '#D55E00',
              'A1003': '#0072B2', 'A1201': '#D55E00', 'A1202': '#33B18F',
              'A1401': '#33B18F', 'A1404': '#DA73B3', 'A1802': '#D55E00',
              'A1803': '#D55E00', 'A1805': '#33B18F', 'A2001': '#D55E00',
              'A2002': '#771853', 'A2003': '#771853', 'A2101': '#33B18F',
              'A2102': '#D55E00', 'A2103': '#0072B2', 'G0101': '#0072B2', 
              'G0102': '#0072B2', 'G0103': '#0072B2', 'G0401': '#0072B2',
              'G0402': '#0072B2', 'G0403': '#33B18F', 'G0404': '#0072B2',
              'G0405': '#0072B2', 'G0406': '#33B18F', 'G0407': '#0072B2',
              'G0408': '#0072B2', 'G0410': '#0072B2', 'G0412': '#0072B2',
              'G0414': '#0072B2', 'G0415': '#0072B2', 'G0417': '#0072B2',
              'G0420': '#33B18F', 'G0601': '#0072B2', 'G0602': '#0072B2',
              'G0702': '#0072B2', 'G0801': '#0072B2', 'G0802': '#0072B2',
              'G0803': '#0072B2', 'G0804': '#0072B2', 'FHON2': '#0072B2',
              'H1B104J': '#0072B2', 'H1B105A': '#0072B2', 'H1B302M': '#0072B2',
              'H2B105J': '#D55E00', 'H3B101A': '#0072B2', 'H3B101J': '#0072B2',
              'H3B101X': '#D55E00', 'H3B102A': '#0072B2', 'H3B102X': '#D55E00',
              'H3B103J': '#0072B2', 'H3B103X': '#D55E00', 'H3B103M': '#0072B2', 
              'H3B104J': '#0072B2', 'H3B104X': '#0072B2', 'H3B107A': '#0072B2', 
              'H3B109M': '#0072B2', 'H3B110M': '#0072B2', 'H3B111A': '#0072B2', 
              'H3B111M': '#0072B2', 'H3B202M': '#D55E00', 'H3B202X': '#0072B2', 
              'H3B203J': '#0072B2', 'H3B203M': '#D55E00', 'H3B204J': '#0072B2', 
              'H3B204M': '#D55E00', 'H3B205J': '#0072B2', 'H3B206M': '#D55E00', 
              'H3B207X': '#0072B2', 'H3B208X': '#0072B2', 'H3B209X': '#33B18F', 
              'H4B101A': '#D55E00', 'H4B102A': '#D55E00', 'H4B103J': '#D55E00', 
              'H4B104A': '#D55E00', 'H4B111J': '#D55E00', 'H4B114J': '#D55E00', 
              'H4B116J': '#D55E00', 'H4B202J': '#0072B2', 'H4B203M': '#D55E00', 
              'H4B204J': '#0072B2', 'H4B205J': '#0072B2', 'H4B206J': '#D55E00', 
              'H4B210M': '#D55E00', 'H4B211M': '#0072B2', 'H4B303J': '#D55E00', 
              'H4B402J': '#0072B2', 'H4B403J': '#D55E00', 'H4B404J': '#D55E00', 
              'H4B405J': '#D55E00', 'H4B406M': '#D55E00', 'H4B410M': '#D55E00', 
              'H4B411M': '#D55E00', 'H4B412M': '#0072B2', 'H4B501J': '#0072B2', 
              'H4B502X': '#0072B2', 'H4B503X': '#0072B2', 'H4B504J': '#33B18F',
              'H4B505J': '#33B18F', 'H4B507J': '#0072B2', 'H4B507X': '#0072B2', 
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': '#D55E00', 
              'DSMZ': '#0072B2', 'HNS-8': 'black', 'GYUN-333': 'black', 
              'FHON13': '#79443B'}

indir = f'{workdir}/trees/RR-TF_HK' #Path to inputs
outdir = f'{workdir}/plots/trees/RR-TF_HK' #Path where outputs will be saved

#Create output directory if it does not exist
if not os.path.exists(outdir):
   os.makedirs(outdir)
   
genes = ['RR', 'HK'] #List of genes for which trees have been generated
suffixes = ['all', 'repset'] #List of datasets (all strains/representative strains)
#Retrieve a list of tree files from gene names and datasets
treefiles = [f'{indir}/{gene}_{suffix}.mafft.faa.treefile' for gene in genes for suffix in suffixes]
   
# =============================================================================
# 2. Loop through tree files and generate output plot
# =============================================================================

for treefile in treefiles:
    outfile = f'{outdir}/{os.path.basename(treefile).split(".")[0]}_support_tree.png' #Set the name of the output
    
    t = Tree(treefile, parser = 'name') #Load the tree file into a tree object
    
    midpoint = t.get_midpoint_outgroup() #Retrieve midpoint
    t.set_outgroup(midpoint) #Use it to root the tree
    
    ts = TreeStyle() #Create tree style object
    ts.show_branch_length = False #Hide branch lengths
    ts.show_branch_support = False #Hide branch supports to add formatted text
    ts.show_leaf_name = False #Hide leaf names to add formatted tex
    ts.scale = 25000 #Set the scale of the tree
    ts.scale_length = 0.01 #Set the length of the legend scale bar
    ts.branch_vertical_margin = 35 #Spacing between branches
    ts.optimal_scale_level = 'full' #Avoid dotted lines to increase the length of branches
    
    ns = NodeStyle() #Create node style
    ns['size'] = 0 #Hide nodes
    ns['vt_line_width'] = 5 #Set width of vertical lines
    ns['hz_line_width'] = 5 #Set width of horizontal lines
    ns['hz_line_type'] = 0 #Horizontal lines will be solid lines
    
    for n in t.traverse(): #Loop through nodes in the tree
       n.set_style(ns) #Apply the style to each node

       if n.name is not None and n not in t.leaves(): #If the node is not a leaf or root node
           support = float(n.name.split('/')[0]) #Retrieve the support value from the node name
           if support >= 50: #If the support value is above 50
               n.support = support #Assign it to the node support property
               if n.support >= 95: #If the node support is above 95
                   color = 'black' #Color the support value in black
               else:
                   color = 'dimgrey' #Otherwise, color it in grey
               support_face = TextFace(int(n.support), fgcolor = color, fsize = 20,
                                       ftype = 'Arial') #Create a text with the support value
               n.add_face(support_face, column = 0, position='branch-top') #Add the text to the corresponding node in the tree
                    
    for leaf in t.leaves(): #Loop through the leaves in the tree
        nleaf = leaf.name #Retrieve the leaf name
        if 'LDX55' in nleaf: #If the leaf name contains this string
            nleaf = nleaf.replace('LDX55', 'IBH001') #Change the string to IBH001
        elif 'APS55' in nleaf: #Same for MP2
            nleaf = nleaf.replace('APS55_RS', 'MP2_')
        elif 'K2W83' in leaf.name: #Same for DSMZ
            nleaf = nleaf.replace('K2W83_RS', 'DSMZ_')
        elif 'MUB42' in leaf.name:
            nleaf = nleaf.replace('MUB42', 'HNS-8')
        elif 'VQ058_RS' in leaf.name:
            nleaf = nleaf.replace('VQ058_RS', 'GYUN-333_')
        else: nleaf = nleaf[3:]
            
        color = leaf_color.get(nleaf.split('_')[0], None) #Retrieve the color getting the strain name from the locus tag
        name_face = TextFace(nleaf, fgcolor = color, ftype = 'Arial', 
                             fsize = 40) #Create a text face with the locus tag, colored by phylogroup 
        leaf.add_face(name_face, column = 0, position = 'branch-right') #Add the locus tag to the leaf
        
    t.ladderize(1) #Change the arrangement of the nodes in the tree so that the root is at the bottom
    t.render(outfile, tree_style = ts, dpi = 400) #, w = 2480, h = 3508) #Save the output plot to PNG
    t.render(outfile.replace('png', 'pdf'), tree_style = ts, dpi = 400) #, w = 2480, h = 3508) #Save the plot to TIFF
    t.render(outfile.replace('png', 'svg'), tree_style = ts, dpi = 400) #, w = 2480, h = 3508) #Save the plot to SVG

end_time = time.time() - start_time #Get total running time of the script
logging.info(f'This script took {end_time/60:2f} minutes.')