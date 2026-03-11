#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021, updated on Tue Mar 3 2026

This script creates plots of the phylogeny of the RR-TF gene and its BLASTp
homologs.
Environment: pixi_phylo/default.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages or modules
# =============================================================================

from ete4 import Tree
from ete4.treeview import TreeStyle, NodeStyle, TextFace, RectFace, CircleFace
import os
import logging, traceback
import sys
import pandas as pd

# =============================================================================
# 0. Logging
# =============================================================================

log = snakemake.log[0] #Log file

with open(log, 'w') as handle: #Open the log file in write mode
    handle.write('') #Overwrite the file

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    '''Function to handle exceptions'''
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger = logging.getLogger() #Create a logger

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ])) #Format to show exceptions

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Define input variables and formatting variables
# =============================================================================

treefile = snakemake.input.tree #Full path to the tree file
outfiles = snakemake.output #Output images with the tree representation
perc_ident = snakemake.input.hits #Full path to the parsed BLASTp results
outdir = os.path.dirname(outfiles[0]) #Output directory

#Dictionary with the colors to use for the different taxonomic groups
color_dict = {'Apilactobacillus': '#D34300', 
            'Lactobacillaceae': '#FA8D33', 
            'Enterococcaceae': '#1B305D', 
            'Carnobacteriaceae': '#9197AE'}

croR = '#F2F1EB' #Color to highlight the genes labelled as CroR

print('Create output directory if it does not exist...')
if not os.path.exists(outdir): #If the output directory does not exist
   os.makedirs(outdir) #Create it
print('Done!')

# =============================================================================
# 1. Load the tree file and set the outgroup
# =============================================================================

print('Load the tree file and root it at its midpoint...')
t = Tree(treefile, parser = 'name') #Load the tree into a Tree object
midpoint = t.get_midpoint_outgroup() #Get midpoint node
t.set_outgroup(midpoint) #Root the tree on the outgroup node
print('Done!')

# =============================================================================
# 2. Create and modify a tree style
# =============================================================================

print('Add formatting to the tree...')
#2.1 Tree formatting
ts = TreeStyle() #Create a tree style 
ts.show_branch_length = False #Hide branch lengths
ts.show_branch_support = False #Hide support values
ts.show_leaf_name = False #Hide leaf names
ts.scale = 2900 #Set the length of the tree
ts.scale_length = 0.1 #Set the length of the tree scale
ts.branch_vertical_margin = 35 #Spacing between branches
ts.optimal_scale_level = 'full' #Avoid dotted lines to increase the length of branches
print('Done!')

print('Add a legend to the tree...')
#2.2 Legend formatting
ts.legend_position = 2 #Legend should be positioned at the upper right
ts.legend.add_face(TextFace('Taxonomy', ftype = 'Arial', fsize = 42), 
                   column = 1) #Add legend title
#Add empty legend next to the title
ts.legend.add_face(TextFace(' ', ftype = 'Arial', fsize = 42), column = 0)

for key in color_dict: #Loop through the keys in the color dictionary
    value = color_dict[key] #Retrieve the color
    ts.legend.add_face(CircleFace(20, value, 'circle'), column =  0) #Create a circle of that color
    ts.legend.add_face(TextFace(f'  {key}', ftype = 'Arial', fsize = 36), 
                       column = 1) #Add the legend text next to the circle

#Add column with the percentage of identity, and an empty column before it
ts.legend.add_face(TextFace('% ID', ftype = 'Arial', fsize = 42), column = 3) 
ts.legend.add_face(TextFace(' ', ftype = 'Arial', fsize = 42), column = 2)

for val in [25, 20, 15]: #Loop through circle sizes
    ts.legend.add_face(CircleFace(val, 'grey', 'circle'), column =  2) #Create a circle of that size
    ts.legend.add_face(TextFace(f'  {val*4}', ftype = 'Arial', fsize = 36), 
                       column = 3) #Add the legend text next to the circle
    
#Add empty legend next to the title again
ts.legend.add_face(TextFace('  ', ftype = 'Arial', fsize = 42), column = 4)

#Add the title of the CroR legend
ts.legend.add_face(TextFace('CroR', ftype = 'Arial', fsize = 42), column = 5)
#Add a rectangle with the color to highlight proteins annotated as CroR
ts.legend.add_face(RectFace(120, 60, croR, croR), column = 5)
print('Done!')

# =============================================================================
# 3. Create and apply a node style
# =============================================================================

print('Create a node style...')
#3.1 Create the node style
ns = NodeStyle() #Create a node style
ns['size'] = 0 #Remove node representations as circles
ns['vt_line_width'] = 1 #Set the width of vertical lines
ns['hz_line_width'] = 1 #Set the width of horizontal lines
ns['hz_line_type'] = 0 #Make horizontal lines solid
ls = ns.copy() #Create another node style to highlight proteins annotated as CroR
ls['bgcolor'] = croR #Assign the background color
print('Done!')

print('Apply the node style and add support values...')
#3.2 Apply the node style on the nodes and add support values
for n in t.traverse(): #Loop through the nodes in the tree
   n.set_style(ns) #Apply the style to the nodes
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
print('Done!')   

# =============================================================================
# 4. Modify the style of the leaf names in the tree
# =============================================================================

print('Add labels and background colors to leaves in the tree')
leaves = t.leaves() #Retrieve the leaves in the tree

perc_id = pd.read_csv(perc_ident, sep = ',') #Dataframe read from the BLASTp results

for leaf in leaves: #Loop through the leaves of the tree
    if 'CroR' in leaf.name: #If the original leaf name includes CroR
        leaf.set_style(ls) #Apply the background highlighting for CroR
    nleaf = leaf.name #Create an additional variable to store the leaf name
    
    prot_id = nleaf.split('.')[0] + '.' + nleaf.split('.')[1][0]
    if prot_id != 'CAI2564285.1':
        id_no = float(perc_id[perc_id['Subject accession'] == prot_id]['% Identity'].iloc[0])
        print(id_no)
    else: id_no = 100

    if '_Apilactobacillus_' in leaf.name: #If Apilactobacillus is in the leaf name
        color = color_dict['Apilactobacillus'] #Assign the corresponding color
    elif ('coccus' in leaf.name and 'Oeno' not in leaf.name and 'Pedio' not in leaf.name): #If it is a species from Enterococcaceae
        color = color_dict['Enterococcaceae'] #Assign the corresponding color
    elif 'Carnobacterium' in leaf.name: #If it is a Carnobacterium
        color = color_dict['Carnobacteriaceae'] #Assign the corresponding color
    else: color = color_dict['Lactobacillaceae'] #Otherwise, use the Lactobacillaceae color
    name_face = TextFace(nleaf.replace('_CroR', '').replace('_', ' ').replace('WP ', 'WP_') + ' ', #Remove CroR from the leaf name and replace underscores with spaces
                         fgcolor = color, fsize = 40, ftype = 'Arial') #Create a text with locus tags
    leaf.add_face(name_face, column = 0, position = 'branch-right') #Add the text to the right leaf in the tree

    #Add the circles with percentages of identity next to the leaves
    perc_face = CircleFace(id_no/4, color, 'circle') #The circles have a size of the % id divided by 4, and are colored by taxonomy
    perc_text = TextFace('100' if id_no == 100 else f'{id_no:.2f}', #Add text face with the % id
                         fgcolor = 'black', fsize = 32, ftype = 'Arial') #Create a text with locus tags
    leaf.add_face(perc_face, column = 1, position = 'branch-right') #Add the circle to the leaf name
    leaf.add_face(perc_text, column = 2, position = 'branch-right') #Add the % id next to the circle
    
print('Done!')

# =============================================================================
# 5. Save the tree to a file
# =============================================================================

print(f'Save the tree to {outfiles[0]} and other image formats...')
#Save the tree to PNG, SVG and PDF
[t.render(outfile, tree_style = ts, dpi = 400) for outfile in outfiles]
print('Done!')