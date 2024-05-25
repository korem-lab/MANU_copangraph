import os
import sys
import tqdm
import math
import networkx as nx
import pandas as pd
import numpy as np
from collections import defaultdict
import Utils.parse_seq as ps
DATA = './data/PanMetagenomeViz/'
RELEVANT_BUGS = [
    's__Enterococcus_B faecium',
    's__Escherichia coli',
    's__Klebsiella pneumoniae'
]
RELEVANT_BUGS_NAME = [e.replace(' ', '_') for e in RELEVANT_BUGS]


def get_seq_length(segments, summary=np.max):
    node_groups = defaultdict(list)
    length = dict()
    for s in segments:
        node_groups[s.nid].append(len(s.seq))
    for k, v in node_groups.items():
        length[k] = summary(v)
    return length


def build_initial_graph(gfa, node_info, boring_weight=0.5):
    
    interesting_node = defaultdict(bool)
    reduction = node_info.loc[:, RELEVANT_BUGS_NAME].apply(lambda x: any(x), axis=1)
    for i in tqdm.tqdm(node_info.index):
        n = node_info.loc[i, 'node']
        interesting_node[n] = interesting_node[n] or reduction[i]
        
    G = nx.DiGraph()
    for e in tqdm.tqdm(gfa):
        if e.type != ps.GFATypes.S:
            continue
        G.add_node(e.nid)
    for e in gfa:
        if e.type != ps.GFATypes.L:
            continue
        G.add_edge(e.l_nid, e.r_nid, weight=1 if (interesting_node[e.l_nid] and interesting_node[e.r_nid]) else boring_weight)
    return G

def expand_graph(G, node_info, node_lengths, keep_self_loops=False, block_size=1000, boring_weight=0.5):
    
    H = G.copy()
    for i in tqdm.tqdm(node_info.index):
        node = node_info.loc[i, 'node']
        if node not in H:
            continue
        # compute record
        label = np.argwhere(node_info.loc[i, RELEVANT_BUGS_NAME].values).flatten()
        if len(label) > 0:
            label = RELEVANT_BUGS_NAME[label[0]]
        else:
            label = '-'
        record = {'label': label, 'num_p': node_info.loc[i, 'num_persistence'], 'node_col': label}
            
        # compute node chain that node will be broken into 
        length = node_lengths[node]
        n_blocks = math.ceil(length/block_size)
        chain = list()
        for i in range(n_blocks):
            chain.append(f'{node}_{i}')
        
        # get adjacent nodes 
        in_adj = list(H.predecessors(node))
        out_adj = list(H.successors(node))
        self_loop = node in out_adj
        
        is_interesting = label != '-'
        w = 1 if is_interesting else boring_weight
        
        # replace node with chain
        H.remove_node(node)
        for n in chain:
            H.add_node(n, **record)
        for i_node in in_adj:
            H.add_edge(i_node, chain[0], weight=w)
        for o_node in out_adj:
            H.add_edge(chain[-1], o_node, weight=w)
        for i in range(len(chain)-1):
            H.add_edge(chain[i], chain[i+1], weight=w)
            
        if keep_self_loops and self_loop:
            G.add_edge(chain[0], chain[1], weight=w)
    return H
        
    
if __name__ == '__main__':
    #G = nx.MultiGraph()
    #G.add_nodes_from([(4, {"Color": "red", "Label": "Klieb"}), (5, {"Color": "green", "Label":"Kleib"})])
    #G.add_edge(4, 5)
    #nx.write_gexf(G, 'testing.gexf')
    
    
    noi = pd.read_pickle(os.path.join(DATA, 'nodes_of_interest.pkl'))
    with open(os.path.join(DATA, 'filtered_gfa.gfa')) as f:
        gfa = list(ps.parse_gfa(f))
    
    G = nx.DiGraph()
    node_lengths = get_seq_length([e for e in gfa if e.type == ps.GFATypes.S])
    
    node_info = noi.loc[noi.node.isin({e.nid for e in gfa if e.type == ps.GFATypes.S})]
    
    print('Build initial graph')
    G = build_initial_graph(gfa, node_info)
    print('Expand graph')
    G = expand_graph(G, node_info, node_lengths)
    print('Write')
    nx.write_gexf(G, os.path.join(DATA, 'filtered_gfa_length_expansion.gexf'))
    
    
#    print(noi.shape)
#    processed_nodes = set()
#    records = list()
#    for i in tqdm.tqdm(filtered_table.index):
#        node = noi.loc[i, 'node']
#        if node not in processed_nodes:
#            node_info = noi.loc[i, :]
#            length = node_lengths[node]
#            label = np.argwhere(node_info[RELEVANT_BUGS_NAME].values).flatten()
#            if len(label) > 0:
#                label = RELEVANT_BUGS_NAME[label[0]]
#            else:
#                label = '-'
#            records += get_records(node, node_info, label, length)
#                
#            
#            record = {
#                "node": node,
#                "num_p": node_info['num_persistence'],
#                "num_c": node_info['num_clearence'],
#                "total": node_info['num_persistence'] + node_info['num_clearence'],
#                "label": label,
#                "node_col": label
#            }
#            records[node] = record
#    G.add_nodes_from(records.items())
#    print('construct edges')
#    for e in tqdm.tqdm(gfa):
#        if e.type == ps.GFATypes.L:
#            G.add_edge(e.l_nid, e.r_nid)
#    for u, v in G.edges:
#        edge_weight = min(G.degree[u], G.degree[v])
#        G[u][v]['weight'] = 1.0/edge_weight
#    nx.write_gexf(G, os.path.join(DATA, 'filtered_gfa_min.gexf'))
    
    
    