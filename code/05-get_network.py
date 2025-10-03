#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:44:38 2023, modified on Thu Oct 2.

Script to create a phylogenetic network.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries and modules
# =============================================================================

from matplotlib import pyplot as plt
import networkx as nx
import os, re
import pandas as pd
import subprocess
import multiprocessing 
import time 
from functools import partial
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline

# =============================================================================
# 1. Function definitions
# =============================================================================

start_time = time.time()

def trim_align(file, log, thr = 1): #Change this to trim the existing MAFFT alignments
    outfile = file.replace('.mafft.faa', '.mafft_trimmed.faa')
    subprocess.run(f'trimal -in {file} -out {outfile} -nogaps -fasta 2> {log}',
                   shell = True)
    return outfile

def needle_align_code(query_seq, target_seq):
    needle_cline = NeedleCommandline(asequence="asis:" + query_seq,
                                      bsequence="asis:" + target_seq,
                                      sprotein=True,
                                      aformat="simple",
                                      gapopen=10,
                                      gapextend=0.5,
                                      outfile='stdout'
                                      )
    out_data, err = needle_cline()
    out_split = out_data.split("\n")
    p = re.compile("\((.*)\)")
    return p.search(out_split[25]).group(1).replace("%", "")

def create_out_tab(out_file: str):
    with open(out_file, 'w') as out:
        out.write('locus1\tlocus2\tPredicted_type\t%identity\n')

def retrieve_subtypes(adh_file: str):
    adh_dict = {}
    with open(adh_file) as tab:
        for line in tab:
            loctags = line.split('\t')
            for i in range(1, len(loctags)):
                locus = loctags[i].strip()
                if 'Adhesin' in locus:
                    continue
                if locus != '-' and i == 1 and 'PLPX' in locus:
                    adh_dict[locus] = 'PLPX'
                elif locus != '-' and i == 1 and 'PLPX' not in locus:
                    adh_dict[locus] = 'CHR1'
                elif locus != '-' and i == 2 and len(locus.split(', ')) == 1:
                    adh_dict[locus] = 'CHR2'
                elif locus != '-' and i == 2 and len(locus.split(', ')) == 2:
                    locus1 = locus.split(', ')[0]
                    locus2 = locus.split(', ')[1]
                    adh_dict[locus1] = 'CHR2'
                    adh_dict[locus2] = 'CHRU'
    return adh_dict
           
def simplify_id(loctag: str):
    new_id = loctag.replace('K2W83_RS', 'DSMZ_').replace('MUB42', 'HNS-8').replace('AKU', '').replace('AAP', '')
    return new_id
             
def process_input(infile: str, out_file: str, adh_dict: str, pos_dict: dict):
    outpath = os.path.dirname(out_file)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    with open(infile) as handle:
        record_list = [record for record in SeqIO.parse(handle, 'fasta') if record.id in adh_dict.keys()]
        for i in range(len(record_list)-1):
            for j in range(i+1, len(record_list)):
                id1 = simplify_id(record_list[i].id)
                id2 = simplify_id(record_list[j].id)
                adh_name = f'{adh_dict[record_list[i].id]}/{adh_dict[record_list[j].id]}'
                perc_ident = needle_align_code(record_list[i].seq, record_list[j].seq)
                pos_dict[(id1, id2)] = perc_ident
                with open(out_file, 'a') as out:
                    out.write(f'{id1}\t{id2}\t{adh_name}\t{perc_ident}\n')

# =============================================================================
# 2. Define the paths to input and output files
# =============================================================================
workdir = os.path.expanduser('~') + '/adhesins' #Working directory
outdir = f'{workdir}/plots/network' #Output directory
adh_file = f'{workdir}/results/adhesin_list.tsv'
log = f'{workdir}/logs/05-network.log'
id_dir = f'{workdir}/results/identity' #File with % of identity

paths = [outdir, id_dir]

[os.makedirs(path) for path in paths if not os.path.exists(path)]

ident = 80 #Threshold percentage of identity
threads = 12

with open(log, 'w') as out:
    out.write('')
    
adh_dict = retrieve_subtypes(adh_file)
adh_list = [f'/{value}' for value in adh_dict.values()]

for suffix in ['', '_repset']:
    aln_file = f'{workdir}/sequences/Muc_adhesins{suffix}.mafft.faa'
    id_file = f'{id_dir}/percentage_identity{suffix}.tab' #Input file with % of identity
    network_out = f'{outdir}/Muc_network{suffix}.svg'
    with open(id_file, 'w') as out:
        out.write('')

    # =============================================================================
    # 3. Run functions to create file with EMBOSS Needle pairwise % of identity
    # =============================================================================

    pos_dict = {}
    aln_out = trim_align(aln_file, log = log)
    process_input(aln_out, id_file, adh_dict, pos_dict)

    # =============================================================================
    # 4. Add colors according to gene category
    # =============================================================================

    color_dict = {'PLPX': '#C6C468', 'CHR1': '#BBE36A', 'CHR2': '#5FB477', 'CHRU': '#8DCC70'} #Dictionary to color nodes by gene type

    color_edge_dict = {'PLPX': '#BEBD8E', 'CHR1': '#BED294', 'CHR2': '#85B193',  
                       'CHRU': '#A4C496'} #Dictionary to color edges by gene type

    with open(id_file) as reader: #Open input file
        df = pd.read_csv(reader, sep = '\t') #Load file contents as a dataframe
        node_names = [] #Empty list to store node information (pairs of nodes)
        edge_list = [] #Empty list to store edge information (pairs of nodes and identity)
        tag_dict = {} #Dictionary to assign a gene type to each node

        for index, row in df.iterrows(): #Loop trough rows in the dataframe
            node_names += [row[0], row[1]] #Add node information
            edge_list.append((row[0], row[1], row[3])) #Add edge information
            tag_dict[row[0]] = row[2].split('/')[0] #Add gene type to first locus tag
            tag_dict[row[1]] = row[2].split('/')[1] #Add gene type to second locus tag

    node_names = set(node_names) #Convert list to set

    G = nx.Graph() #Create the network graph
    G.add_nodes_from(node_names) #Add nodes to graph
    G.add_weighted_edges_from(edge_list) #Add edges to graph, considering weights (% identity)

    colors = [color_dict[tag_dict[loctag]] for loctag in node_names] #Get the colors of each node

    edge_colors = ['#d1d1d1' if tag_dict[edge[0]] != tag_dict[edge[1]] else color_edge_dict[tag_dict[edge[0]]] for edge in list(G.edges())] #Get edge colors (grey if the genes belong to different subtypes)

    # =============================================================================
    # 3.Plot everything as a network (using the code above)
    # =============================================================================

    fig, ax = plt.subplots() #Create plot
    ax.margins(0.1) #Set plot margins
    ax.axis('off') #Remove the plot axis
    fig.set_size_inches(48, 32) #Set figure size

    weights = [pos_w[2]['weight']/10 if pos_w[2]['weight'] > ident else 0 for pos_w in list(G.edges.data())] #Adjust the weights of the edges and remove edges with identity < 70%
    distances = dict(nx.shortest_path_length(G, weight='weight')) #Convert the weights to distances
    for k, v in distances.items(): #Loop through node names (k = node 1, v = node 2 and distance)
        for k2, v2 in v.items(): #Loop through node names and distances (k2 = node 2, v2 = distance)
            print(k, k2) #Print the name of both nodes
            if distances[k][k2] != 0: #If the distance is not 0
                distances[k][k2] = 1/distances[k][k2] #Reverse the distance (so that the higher the % identity, the lower the distance)

    net_pos = nx.kamada_kawai_layout(G, dist = distances) #Generate network

    # Divide this into a command to plot edges with lower alpha and a command to plot
    # nodes with higher alpha

    nx.draw_networkx_edges(G, pos = net_pos, edgelist = G.edges(), alpha = 0.5,
                            edge_color = edge_colors, width = weights, ax = ax) #Plot edges

    nx.draw_networkx_nodes(G, pos = net_pos, nodelist = G.nodes(), alpha = 1, 
                           node_color = colors, node_shape = '*', 
                           node_size = 1200, ax = ax) #Plot nodes

    nx.draw_networkx_labels(G, pos=net_pos, font_size = 20, font_weight='bold',
                            verticalalignment = 'baseline', ax = ax,
                            horizontalalignment = 'center') #Add labels to nodes


    # =============================================================================
    # 4. Here I save the network to files
    # =============================================================================

    fig.savefig(network_out, format='svg', dpi=800, pad_inches = 0) #Save network to SVG
    fig.savefig(network_out.replace('.svg', '.png'), format='png', dpi=800, 
                pad_inches = 0) #Save network to PNG
    fig.savefig(network_out.replace('.svg', '.tiff'), format='tiff', dpi=800, 
                pad_inches = 0) #Save netwoek to TIFF
    fig.savefig(network_out.replace('.svg', '.pdf'), format='pdf', dpi=800, 
                pad_inches = 0) #Save netwoek to PDF

    end_time = time.time() - start_time
    print(f'This script took {end_time/60:2f} minutes.')