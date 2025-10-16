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
import time
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
import multiprocessing
from functools import partial

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
    needle_cline = NeedleCommandline(asequence='asis:' + query_seq,
                                      bsequence='asis:' + target_seq,
                                      sprotein=True,
                                      aformat='simple',
                                      gapopen=10,
                                      gapextend=0.5,
                                      outfile='stdout'
                                      )
    out_data, err = needle_cline()
    out_split = out_data.split('\n')
    p = re.compile('\((.*)\)')
    return p.search(out_split[25]).group(1).replace('%', '')

def create_out_tab(out_file: str):
    with open(out_file, 'w') as out:
        out.write('locus1\tlocus2\tPredicted_type\t%identity\n')

def retrieve_subtypes(adh_file: str, adhesin: str):
    adh_dict = {}
    with open(adh_file) as tab:
        for line in tab:
            loctags = line.split('\t')
            if len(loctags) == 2:
                locus = loctags[1].strip()
                adh_dict[locus] = adhesin
            else:
                for i in range(1, len(loctags)):
                    locus = loctags[i].strip()
                    if 'Adhesin' in locus:
                        continue
                    if adhesin == 'MucBP+LPXTG':
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
                    elif adhesin == 'MubB2+LPXTG' and locus != '-':
                        adh_dict[locus] = f'MUB{i}'

    return adh_dict
           
def simplify_id(loctag: str):
    new_id = loctag.replace('K2W83_RS', 'DSMZ_').replace('MUB42', 'HNS-8').replace('LDX55', 'IBH001').replace('APS55_RS', 'MP2_').replace('VQ058_RS', 'GYUN-333_').replace('AKU', '').replace('AAP', '')
    return new_id
             
def get_length_dict(faa_file): #Function not used
    len_dict = {}
    with open(faa_file) as faa:
        for record in SeqIO.parse(faa, 'fasta'):
            len_dict[record.id] = len(record.seq)
    return len_dict

def perc_ident(aln_file, len_dict, adh_dict, outfile):
    
    outpath = os.path.dirname(outfile)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    per_ident = {}
    with open(aln_file) as aln:
        aln_read = list(SeqIO.parse(aln, 'fasta'))
        for i in range(0, len(aln_read)-1):
            record1 = aln_read[i]
            for j in range(1, len(aln_read)):
                record2 = aln_read[j]
                id_count = 0
                min_len = min(len_dict[record1.id], len_dict[record2.id])
                for k in range(0, len(record1.seq)):
                    if not (record1.seq[k] == '-' and record2.seq[k] == '-') and record1.seq[k] == record2.seq[k]:
                        id_count += 1
                        
                per_ident[(record1.id, record2.id)] = (id_count/min_len)*100
                
    with open(outfile, 'w') as out_handle:
        out_handle.write('locus1\tlocus2\tsubtypes\t%id\n')
        [out_handle.write(f'{simplify_id(locid[0])}\t{simplify_id(locid[1])}\t{adh_dict[locid[0]]}/{adh_dict[locid[1]]}\t{per_ident[(locid[0], locid[1])]}\n') for locid in per_ident.keys()]
    return per_ident

def full_process(adhesin_suffix: tuple, adh_path: str, id_dir: str, adh_list_dir: str, log: str):
    adhesin = adhesin_suffix[0]
    suffix = adhesin_suffix[1]
    in_aln = f'{adh_path}/{adhesin}{suffix}.mafft.faa'
    infile = trim_align(in_aln, log = log)
    out_file = f'{id_dir}/{adhesin}_percentage_identity{suffix}.tab'
    adh_file = f'{adh_list_dir}/{adhesin}_list.tsv'
    adh_dict = retrieve_subtypes(adh_file, adhesin)
    
    outpath = os.path.dirname(out_file)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    create_out_tab(out_file)
        
    with open(infile) as handle:
        record_list = [record for record in SeqIO.parse(handle, 'fasta') if record.id in adh_dict.keys()]
        for i in range(len(record_list)-1):
            for j in range(i+1, len(record_list)):
                id1 = simplify_id(record_list[i].id)
                id2 = simplify_id(record_list[j].id)
                adh_name = f'{adh_dict[record_list[i].id]}/{adh_dict[record_list[j].id]}'
                perc_ident = needle_align_code(record_list[i].seq, record_list[j].seq)
                # pos_dict[(id1, id2)] = perc_ident
                with open(out_file, 'a') as out:
                    out.write(f'{id1}\t{id2}\t{adh_name}\t{perc_ident}\n')
    return out_file

def create_network(id_file):
    adhesin_type = os.path.basename(id_file).split('_')[0]
    if 'repset' in id_file:
        suff = '_repset'
    else: suff = ''
    network_out = f'{outdir}/{adhesin_type}{suff}.svg'
    with open(id_file) as reader: #Open input file
    
        #Segment to create the network
        
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

       #Segment to plot the network

        fig, ax = plt.subplots() #Create plot
        ax.margins(0.1) #Set plot margins
        ax.axis('off') #Remove the plot axis
        fig.set_size_inches(48, 32) #Set figure size

        weights = [pos_w[2]['weight']/10 if pos_w[2]['weight'] > ident else 0 for pos_w in list(G.edges.data())] #Adjust the weights of the edges and remove edges with identity < 70%
        distances = dict(nx.shortest_path_length(G, weight='weight')) #Convert the weights to distances
        for k, v in distances.items(): #Loop through node names (k = node 1, v = node 2 and distance)
            for k2, v2 in v.items(): #Loop through node names and distances (k2 = node 2, v2 = distance)
                # print(k, k2) #Print the name of both nodes
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

        fig.savefig(network_out, format = 'svg', dpi = 800, pad_inches = 0) #Save network to SVG
        fig.savefig(network_out.replace('.svg', '.png'), format = 'png', dpi = 800, 
                    pad_inches = 0) #Save network to PNG
        fig.savefig(network_out.replace('.svg', '.tiff'), format = 'tiff', dpi = 800, 
                    pad_inches = 0) #Save network to TIFF
        fig.savefig(network_out.replace('.svg', '.pdf'), format = 'pdf', dpi = 800, 
                    pad_inches = 0) #Save network to PDF
        
    return network_out
                    
color_dict = {'PLPX': '#C6C468', 'CHR1': '#BBE36A', 'CHR2': '#5FB477', 
              'CHRU': '#8DCC70', 'SH3b': '#037971', 'Gtf2': '#5BC0EB', 
              'collagen-binding': '#3A1772', 'MUB1': '#DF2935', 
              'MUB2': '#FF9D3B', 'MUB3': '#FDCA40'} #Dictionary to color nodes by gene type

color_edge_dict = {'PLPX': '#BEBD8E', 'CHR1': '#BED294', 'CHR2': '#85B193',  
                   'CHRU': '#A4C496', 'SH3b': '#51A9A3', 'Gtf2': '#9DCBDF', 
                   'collagen-binding': '#6F46B0', 'MUB1': '#D27E84', 
                   'MUB2': '#E7BA8D', 'MUB3': '#E6CF8E'} #Dictionary to color edges by gene type

# =============================================================================
# 2. Define the paths to input and output files
# =============================================================================

workdir = os.path.expanduser('~') + '/adhesins' #Working directory
outdir = f'{workdir}/plots/network' #Output directory
log = f'{workdir}/logs/05-network.log'
id_dir = f'{workdir}/results/identity' #File with % of identity
adh_path = f'{workdir}/alignments/adhesins'
adh_list_dir = f'{workdir}/results/adhesin_lists'

paths = [outdir, id_dir, adh_path]

[os.makedirs(path) for path in paths if not os.path.exists(path)]

ident = 80 #Threshold percentage of identity
threads = 48
colors = ['green', 'blue', 'red', 'gray', 'orange']
adh_types = ['MucBP+LPXTG', 'MubB2+LPXTG', 'Gtf2', 'collagen-binding', 'SH3b']

with open(log, 'w') as out:
    out.write('')
    
# =============================================================================
# 3. Run functions to create file with EMBOSS Needle pairwise % of identity
# with multithreading
# =============================================================================

adh_suffix = [(adh, suffix) for adh in adh_types for suffix in ['', '_repset']]

if __name__ == '__main__': 
    pool = multiprocessing.Pool() 
    pool = multiprocessing.Pool(processes=threads)
    outputs = pool.map(partial(full_process, adh_path = adh_path, id_dir = id_dir, adh_list_dir = adh_list_dir, log = log), adh_suffix)
    print("Input: {}".format(adh_suffix))
    print("Output: {}".format(outputs))
    
end_time = time.time() - start_time
print(f'So far, the script has taken {(end_time/60):2f} minutes.')
    
# =============================================================================
# 4. Generate the network plots with multithreading
# =============================================================================

if __name__ == '__main__': 
    pool = multiprocessing.Pool() 
    pool = multiprocessing.Pool(processes=threads)
    networks = pool.map(create_network, outputs)
    print("Input: {}".format(outputs))
    print("Output: {}".format(networks))

end_time = time.time() - start_time
print(f'So far, the script has taken {(end_time/60):2f} minutes.')