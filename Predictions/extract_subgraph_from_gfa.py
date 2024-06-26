import parse_seq 
from scipy import sparse as sp
import sys 
import pandas as pd
import os

def build_adjM(nodes, edges):

    nodes = sorted(list(set(nodes)))
    node_to_index = {v: i for i, v in enumerate(nodes)}
    num_nodes = len(nodes)
    adjM = sp.lil_matrix((num_nodes, num_nodes), dtype=bool)
    set_nodes = set(nodes) 
    for u, v in edges:
        if u not in set_nodes or v not in set_nodes:
            continue
        adjM[node_to_index[u], node_to_index[v]] = True
        adjM[node_to_index[v], node_to_index[u]] = True

    adjM_df = pd.DataFrame.sparse.from_spmatrix(adjM, index=nodes, columns=nodes)
    return adjM_df

def gfa2adjM(fl):

    nodes = [e for e in parse_seq.parse_gfa(open(fl)) if e.type == parse_seq.GFATypes.S]
    edges = [e for e in parse_seq.parse_gfa(open(fl)) if e.type == parse_seq.GFATypes.L]
    edges = [(e.l_nid, e.r_nid) for e in edges]
    node_ids = set([e.nid for e in nodes])
    contigs = [(n.nid, n.seq) for n in nodes]
    return build_adjM(node_ids, edges), contigs

def dfs(adjM, nodes, u, depth):
    
    nodes += {u}
    if depth == 0:
        return
    
    adj_nodes = set(adjM.index[adjM[u,:]])
    nodes += adj_nodes
    for v in adj_nodes:
        dfs(adjM, nodes, v, depth-1)
        
def write_gfa_subgraph(adjM, nodes_lookup, subgraph_nodes, gfa_name):
    with open(gfa_name, 'w') as f:
        for n in subgraph_nodes:
            f.write(nodes_lookup[n] + '\n')
        for u in subgraph_nodes:
            for v in subgraph_nodes:
                if adjM[u, v]:
                    f.write(f'{u}\t+\t{v}\t+\t0M')
                
if __name__ == '__main__':

    if len(sys.argv) != 5:
        print('Usage: <exe> <gfa> <node> <distance> <out_dir>')
        sys.exit()
    
    links = list()
    adjM = gfa2adjM(sys.argv[1])
    
    start_node = sys.argv[2]
    distance = sys.argv[3]
    assert start_node in adjM.index
    
    nodes_lookup = {e.name: e for e in parse_seq.parse_gfa(open(sys.argv[1])) if e.type == parse_seq.GFATypes.S}
    for d in range(distance):
        subgraph_nodes = set()
        dfs(adjM, subgraph_nodes, start_node, distance)
        write_gfa_subgraph(nodes_lookup, subgraph_nodes, f'{start_node}_{d}.gfa')
    
        
    
        