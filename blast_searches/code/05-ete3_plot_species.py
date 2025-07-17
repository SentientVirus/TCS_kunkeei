#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates plots of the domains of GH70 and GH32 proteins
next to a strain phylogeny. To plot a specific gene type, the script has
to be manually updated.
Environment: ete3.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages or modules
# =============================================================================
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import os

# =============================================================================
# 1. Define input variables and formatting variables
# =============================================================================

# #Dictionary to assign colors to each locus tag depending on strain phylogroup
# leaf_color = {'A0901': '#D55E00', 'A1001': '#771853', 'A1002': '#D55E00',
#               'A1003': '#0072B2', 'A1201': '#D55E00', 'A1202': '#33B18F',
#               'A1401': '#33B18F', 'A1404': '#DA73B3', 'A1802': '#D55E00',
#               'A1803': '#D55E00', 'A1805': '#33B18F', 'A2001': '#D55E00',
#               'A2002': '#771853', 'A2003': '#771853', 'A2101': '#33B18F',
#               'A2102': '#D55E00', 'A2103': '#0072B2', 'G0101': '#0072B2', 
#               'G0102': '#0072B2', 'G0103': '#0072B2', 'G0401': '#0072B2',
#               'G0402': '#0072B2', 'G0403': '#33B18F', 'G0404': '#0072B2',
#               'G0405': '#0072B2', 'G0406': '#33B18F', 'G0407': '#0072B2',
#               'G0408': '#0072B2', 'G0410': '#0072B2', 'G0412': '#0072B2',
#               'G0414': '#0072B2', 'G0415': '#0072B2', 'G0417': '#0072B2',
#               'G0420': '#33B18F', 'G0601': '#0072B2', 'G0602': '#0072B2',
#               'G0702': '#0072B2', 'G0801': '#0072B2', 'G0802': '#0072B2',
#               'G0803': '#0072B2', 'G0804': '#0072B2', 'FHON2': '#0072B2',
#               'H1B104J': '#0072B2', 'H1B105A': '#0072B2', 'H1B302M': '#0072B2',
#               'H2B105J': '#D55E00', 'H3B101A': '#0072B2', 'H3B101J': '#0072B2',
#               'H3B101X': '#D55E00', 'H3B102A': '#0072B2', 'H3B102X': '#D55E00',
#               'H3B103J': '#0072B2', 'H3B103X': '#D55E00','H3B103M': '#0072B2', 
#               'H3B104J': '#0072B2', 'H3B104X': '#0072B2', 'H3B107A': '#0072B2', 
#               'H3B109M': '#0072B2', 'H3B110M': '#0072B2', 'H3B111A': '#0072B2', 
#               'H3B111M': '#0072B2', 'H3B202M': '#D55E00', 'H3B202X': '#0072B2', 
#               'H3B203J': '#0072B2', 'H3B203M': '#D55E00', 'H3B204J': '#0072B2', 
#               'H3B204M': '#D55E00', 'H3B205J': '#0072B2', 'H3B206M': '#D55E00', 
#               'H3B207X': '#0072B2', 'H3B208X': '#0072B2', 'H3B209X': '#33B18F', 
#               'H4B101A': '#D55E00', 'H4B102A': '#D55E00', 'H4B103J': '#D55E00', 
#               'H4B104A': '#D55E00', 'H4B111J': '#D55E00', 'H4B114J': '#D55E00', 
#               'H4B116J': '#D55E00', 'H4B202J': '#0072B2', 'H4B203M': '#D55E00', 
#               'H4B204J': '#0072B2', 'H4B205J': '#0072B2', 'H4B206J': '#D55E00', 
#               'H4B210M': '#D55E00', 'H4B211M': '#0072B2', 'H4B303J': '#D55E00', 
#               'H4B402J': '#0072B2', 'H4B403J': '#D55E00', 'H4B404J': '#D55E00', 
#               'H4B405J': '#D55E00', 'H4B406M': '#D55E00', 'H4B410M': '#D55E00', 
#               'H4B411M': '#D55E00', 'H4B412M': '#0072B2', 'H4B501J': '#0072B2', 
#               'H4B502X': '#0072B2', 'H4B503X': '#0072B2', 'H4B504J': '#33B18F',
#               'H4B505J': '#33B18F', 'H4B507J': '#0072B2', 'H4B507X': '#0072B2', 
#               'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': '#D55E00', 
#               'DSM': '#0072B2'}

treefile = os.path.expanduser('~') + '/mucoid_project/blast_searches/fasta/blastp_formatted.mafft.faa.treefile' #Full path to the tree file
outfile = os.path.expanduser('~') + '/mucoid_project/blast_searches/plots/blastp_multispecies_tree_support.png' #Output image with the tree representation
outdir = os.path.dirname(outfile) #Output directory

if not os.path.exists(outdir): #If the output directory does not exist
   os.makedirs(outdir) #Create it

# =============================================================================
# 1. Load the tree file and set the outgroup
# =============================================================================

t = Tree(treefile, format = 0) #Load the tree into a Tree object
midpoint = t.get_midpoint_outgroup()
#outnode = t.get_common_ancestor('AOR73699.1_4_3_GTFB_Limosilactobacillus_fermentum', 'AAU08014.2_4_6_GTFB_Limosilactobacillus_reuteri') #Retrieve the outgroup node
t.set_outgroup(midpoint) #Root the tree on the outgroup node

# =============================================================================
# 2. Create and modify a tree style
# =============================================================================

ts = TreeStyle() #Create a tree style 
ts.show_branch_length = False #Hide branch lengths
ts.show_branch_support = False #Hide support values
ts.show_leaf_name = False #Hide leaf names
ts.scale = 2000 #Set the scale of the tree
ts.scale_length = 0.2 #Set the length of the tree
# ts.tree_width = 1 #Set the width of the branches, overwritten by the scale
# ts.force_topology = True #Forces all the branches to be of equal length

# =============================================================================
# 3. Create and apply a node style
# =============================================================================

ns = NodeStyle() #Create a node style
ns['size'] = 0 #Remove node representations as circles
ns['vt_line_width'] = 1 #Set the width of vertical lines
ns['hz_line_width'] = 1 #Set the width of horizontal lines
ns['hz_line_type'] = 0 #Make horizontal lines solid

for n in t.traverse(): #Loop through the nodes in the tree
   n.set_style(ns) #Apply the style to the nodes
   if n not in t.get_leaves() and n.support >= 50: #If the node is not a leaf and the support value is higher of equal to 50
       if n.support >= 95: #If the support value is > 95
           color = 'black' #Color the support value in black
       else: color = 'dimgrey' #Otherwise, color it in grey

       support_face = TextFace(int(n.support), fgcolor = color, fsize = 24,
                               ftype = 'Arial') #Create a text with the support value
       n.add_face(support_face, column = 0, position='branch-top') #Add the text to the corresponding node in the tree
   
# =============================================================================
# 4. Modify the style of the leaf names in the tree
# =============================================================================

leaves = t.get_leaves() #Sort the leaves (locus tags) by phylogeny

for leaf in leaves: #Loop through the leaves of the tree
    nleaf = leaf.name #Create an additional variable to store the leaf name

    if '_Apilactobacillus_' in leaf.name:
        color = '#F1690E'
    # elif 'Limosilactobacillus_reuteri' in leaf.name or 'L_fermentum' in leaf.name:
    #     color = '#350CA6'
    # elif 'Leuconostoc_mesenteroides' in leaf.name or 'L_citreum' in leaf.name:
    #     color = '#7F0797'
    elif '_Streptococcus_' in leaf.name:
        color = '#09569A'
    elif ('coccus' in leaf.name and 'Oeno' not in leaf.name and 'Pedio' not in leaf.name) or 'Erwinia' in leaf.name:
        color = '#0B7551'
    # elif '_P_myokoensis' in leaf.name or '_F_' in leaf.name or '_Nicoliella_' in leaf.name or '_Convinina_' in leaf.name:
    #     color = '#52312E'
    else: color = '#CF0D61'
    name_face = TextFace(nleaf, fgcolor = color, fsize = 40, ftype = 'Arial') #Create a text with locus tags
    leaf.add_face(name_face, column = 0, position = 'branch-right') #Add the text to the right leaf in the tree
    
# pairs = {'GS1': ('H4B505J_12880', 'A1003_12540'), 
#          'GS2': ('H4B504J_13480', 'H4B505J_12900'),
#          'GS3': ('H4B406M_13450', 'H4B111J_13560'),
#          'BRS': ('H4B505J_12890', 'LDX55_06330'),
#          'NGB': ('H4B504J_05510', 'WP_353317313.1_Apilactobacillus_apinorum')}
# # colors = ['#336E74', '#FF707C', '#FFE570', '#FF986F', '#FFC36F', '#BDE384', 
# #           '#84E8A7', '#656ED4', '#C870FF', '#F87BFF', '#84CDE8']
# colors = ['#FF707C', '#FFE570', '#FF986F', '#84E8A7', '#336E74'] #, '#FFC36F']
# node_list = []
# node_content = t.get_cached_content()
# for key in pairs.keys():
#     to_collapse = t.get_common_ancestor(pairs[key][0], pairs[key][1])
#     node_list.append(to_collapse)
#     to_collapse.name = key
    
# for i in range(len(node_list)):
#     ns2 = NodeStyle() #Create a node style
#     ns2['shape'] = 'sphere'
#     ns2['size'] = len(node_content[node_list[i]]) * 10
#     ns2['fgcolor'] = colors[i]
#     ns2['draw_descendants'] = False
#     node_list[i].set_style(ns2)
    
# GS_text = TextFace('GS', fsize = 36, ftype = 'Arial')
# BRS_text = TextFace('BRS', fsize = 36, ftype = 'Arial')
# texts = [GS_text, BRS_text]
# text_colors = ['#F7CBCB', '#D0F7CB'] #'#D2CCF8']
# for j in range(len(texts)):
#     texts[j].margin_top = 0
#     texts[j].margin_right = 20
#     texts[j].margin_left = 20
#     texts[j].margin_bottom = 0
#     texts[j].background.color = text_colors[j]
    
# GS_nodes = [node for node in t.get_leaves() if 'SR_' in node.name]
# BRS_nodes = [node for node in t.get_leaves() if 'BRS_' in node.name]
# for n in GS_nodes:
#     n.add_face(GS_text, column = 0, position = 'branch-top')
    
# for n in BRS_nodes:
#     n.add_face(BRS_text, column = 0, position = 'branch-top')

    
# =============================================================================
# 5. Change final details of tree formatting and save the tree to a file
# =============================================================================

t.ladderize(1) #Reverse the order in which the leaves of the tree are plotted
# GS_BRS = t.get_common_ancestor('H1B302M_12890', 'WP_249510840.1_Apilactobacillus_apisilvae')
# GS_BRS.ladderize(0)
# GS1 = t.get_common_ancestor('WP_353322098.1_Apilactobacillus_apinorum', 'WP_260116903.1_Nicoliella_spurrieriana')
# GS1.ladderize(1)
# GSA = t.get_common_ancestor('BAA14241.1_MSR_Streptococcus_sobrinus', 'WP_323728076.1_Streptococcus_sp_G418')
# GSA.ladderize(1)
# GSB = t.get_common_ancestor('WP_186432386.1_Oenococcus_oeni', 'WP_261912693.1_Lentilactobacillus_hilgardii')
# GSB.ladderize(1)
# GSC = t.get_common_ancestor('WP_273750349.1_Leuconostoc_mesenteroides', 'WP_172692410.1_Latilactobacillus_sakei')
# GSC.ladderize(1)
# GSD = t.get_common_ancestor('WP_010690266.1_Ligilactobacillus_animalis', 'WP_273774597.1_Limosilactobacillus_reuteri')
# GSD.ladderize(1)
# BRS = t.get_common_ancestor('WP_080732772.1_Leuconostoc_mesenteroides', 'WP_248720543.1_Convivina_intestini')
# BRS.ladderize(1)
# t.convert_to_ultrametric() #Convert the tree to a cladogram
t.render(outfile, tree_style = ts) #Save the tree to a PNG image
t.render(outfile.replace('png', 'svg'), tree_style = ts) #Save the tree to an SVG image
t.render(outfile.replace('png', 'pdf'), tree_style = ts) #Save the tree to a PDF image
