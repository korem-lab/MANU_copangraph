import networkx as nx
import Utils.parse_seq as ps
import pandas as pd
import numpy as np
import os
from collections import defaultdict
import tqdm
import sys
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


def build_initial_graph(gfa):
    G = nx.Digraph()
    
    for e in gfa:
        if e.type != ps.GFATypes.S:
            continue
        G.add_node(e.nid)
    for e in gfa:
        if e.type != ps.GFATypes.L:
            continue
        G.add_edge(e.l_nid, e.r_nid)
    return G

def expand(G, node, record):
    in_adj = list(G.predecessors(node))
    out_adj = list(G.successors(node))
    n_blocks = 10
    chain = list()
    for i in range(n_blocks):
        chain.append(f'{node}_{i}')
    G.remove_node(node)
    for n in chain:
        G.add_node(n, **record)
    for i_node in in_adj:
        G.add_edge(i_node, chain[0])
    for o_node in out_adj:
        G.add_edge(chain[-1], o_node)
    for i in range(len(chain)-1):
        G.add_edge(chain[i], chain[i+1])
    return G

def expand_graph(G, node_info, node_lengths, keep_self_loops=False, block_size=1000):
    
    H = G.copy()
    for i in tqdm.tqdm(node_info.index):
        node = node_info.loc[i, 'node']
        
        # compute record
        label = np.argwhere(node_info[RELEVANT_BUGS_NAME].values).flatten()
        if len(label) > 0:
            label = RELEVANT_BUGS_NAME[label[0]]
        else:
            label = '-'
        record = {'label': label, 'num_p': node_info['num_persistence']}
            
        # compute node chain that node will be broken into 
        length = node_lengths[node]
        n_blocks = math.ceil(length/block_size)
        chain = list()
        for i in range(n_blocks)
            chain.append(f'{node}_{i}')
        
        # get adjacent nodes 
        in_adj = list(H.predecessors(node))
        out_adj = list(H.successors(node))
        self_loop = node in out_adj
        
        # replace node with chain
        for n in chain:
            H.add_node(n, **record)
        for i_node in in_adj:
            G.add_edge(i_node, chain[0])
        for o_node in out_adj:
            G.add_edge(chain[-1], o_node)
        for i in range(len(chain)-1):
            G.add_edge(chain[i], chain[i+1])
            
        if keep_self_loops and self_loop:
            G.add_edge(chain[0], chain[1])
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
    
    print('Build initial graph')
    node_info = noi.loc[noi.node.isin({e.nid for e in gfa if e.type == ps.GFATypes.S})]
    
    G = build_initial_graph(gfa)
    G = expand_graph(G, node_info, node_lengths)
    
    # filter down to relevant nodes
    
    
    print(noi.shape)
    processed_nodes = set()
    records = list()
    for i in tqdm.tqdm(filtered_table.index):
        node = noi.loc[i, 'node']
        if node not in processed_nodes:
            node_info = noi.loc[i, :]
            length = node_lengths[node]
            label = np.argwhere(node_info[RELEVANT_BUGS_NAME].values).flatten()
            if len(label) > 0:
                label = RELEVANT_BUGS_NAME[label[0]]
            else:
                label = '-'
            records += get_records(node, node_info, label, length)
                
            
            record = {
                "node": node,
                "num_p": node_info['num_persistence'],
                "num_c": node_info['num_clearence'],
                "total": node_info['num_persistence'] + node_info['num_clearence'],
                "label": label,
                "node_col": label
            }
            records[node] = record
    G.add_nodes_from(records.items())
    print('construct edges')
    for e in tqdm.tqdm(gfa):
        if e.type == ps.GFATypes.L:
            G.add_edge(e.l_nid, e.r_nid)
    for u, v in G.edges:
        edge_weight = min(G.degree[u], G.degree[v])
        G[u][v]['weight'] = 1.0/edge_weight
    nx.write_gexf(G, os.path.join(DATA, 'filtered_gfa_min.gexf'))
    
    
    