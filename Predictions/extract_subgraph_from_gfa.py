import sys 
import os
import re
import pandas as pd
import parse_seq 
from scipy import sparse as sp
from collections import defaultdict
COPAN_DIR = '/manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/analyses/mdro_prediction/scratch/copangraphs/'
OUT_DIR = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/mdro'

def edge_kv(u, v, uori, vori):
    if u < v:
        return (u, v, uori, vori)
    else:
        return (v, u, (vori+1)%2, (uori+1)%2)
    
def build_adjM(nodes, edges):

    nodes = sorted(list(set(nodes)))
    set_nodes = set(nodes) 
    adjM = defaultdict(set)
    edge_lookup = defaultdict(set) 
    for u, v, uori, vori in edges:
        if u not in set_nodes or v not in set_nodes:
            continue
        u, v, uori, vori = edge_kv(u, v, uori, vori)
        edge_lookup[(u, v)].add((u, v, uori, vori))
        adjM[u].add(v)
        adjM[v].add(u)
    return adjM, edge_lookup

def gfa2adjM(fl):
    nodes = [e for e in parse_seq.parse_gfa(open(fl)) if e.type == parse_seq.GFATypes.S]
    edges = [e for e in parse_seq.parse_gfa(open(fl)) if e.type == parse_seq.GFATypes.L]
    edges = [(e.l_nid, e.r_nid, e.l_ori, e.r_ori) for e in edges]
    node_ids = set([e.nid for e in nodes])
    return build_adjM(node_ids, edges)

def dfs(adjM, nodes, u, halt_degree, depth=5, halt=False):
    u = str(u)
    nodes.update({u})
    if depth == 0 or halt:
        return
    
    #adj_nodes = set(adjM.index[adjM.loc[str(u),:]])
    adj_nodes = adjM[u]
    halt = len(adj_nodes) >= halt_degree
    nodes.update(adj_nodes)
    for v in adj_nodes:
        dfs(adjM, nodes, v, halt_degree, depth-1, halt)
        
def write_gfa_subgraph(adjM, nodes_lookup, edge_lookup, subgraph_nodes, gfa_name):
    with open(gfa_name, 'w') as f:
        for n in subgraph_nodes:
            for e in nodes_lookup[n]:
                e.write(f)
        for u in subgraph_nodes:
            for v in subgraph_nodes:
                if v in adjM[u] or u in adjM[v]:
                    if u < v:
                        for _u, _v, uori, vori in edge_lookup[(u, v)]:
                            uori = '+' if uori==1 else '-'
                            vori = '+' if vori==1 else '-'
                            f.write(f'{_u}\t{uori}\t{_v}\t{vori}\t0M\n')
                    else:
                        for _u, _v, uori, vori in edge_lookup[(u, v)]:
                            uori = '+' if uori==1 else '-'
                            vori = '+' if vori==1 else '-'
                            f.write(f'{_u}\t{uori}\t{_v}\t{vori}\t0M\n')

def write_contigs(sample_names, nodes_lookup, subgraph_nodes, fname):
    
    contigs = list()
    get_sample = lambda x: int(x.rest[0].split(':')[3])
    get_cn = lambda x: x.rest[0].split(':')[2]
    
    # get segments grouped by sample
    segments = defaultdict(list)
    for n in subgraph_nodes:
        for s in nodes_lookup[n]:
            segments[get_sample(s)].append(s)
            
    print('num_samples: ', len(segments)) 
    for k, v in segments.items():
        sample_name = sample_names.loc[k, 0]
        cset = {get_cn(e) for e in v}
        with open(os.path.join(OUT_DIR, f'{sample_name}.final.contigs.fa')) as f:
            for fa in parse_seq.parse(f, parse_seq.Fasta):
                cn = fa.hdr.split(' ')[0]
                if cn in cset:
                    fa.hdr = f'{cn} {k} {sample_name}'
                    contigs.append(fa)
    print(len(contigs))
    with open(fname, 'w') as f:
        for c in contigs:
            c.write(f)
        
    
if __name__ == '__main__':

    if len(sys.argv) != 6:
        print('Usage: <exe> <copan_prefix> <top_features> <halting_degree> <out_dir> <contigs>')
        sys.exit()
    
    # load gfa 
    links = list()
    copan_pref = sys.argv[1]
    copangraph = os.path.join(COPAN_DIR, copan_pref + '.gfa')
    emap = os.path.join(COPAN_DIR, copan_pref + '.ecolor.feature_map.csv')
    nmap = os.path.join(COPAN_DIR, copan_pref + '.ncolor.feature_map.csv')
    adjM, edge_lookup = gfa2adjM(copangraph)
    sample_names = pd.read_csv(sys.argv[5], header=None)
    
    # load maps
    nmap = pd.read_csv(nmap, header=None)
    emap = pd.read_csv(emap, header=None)
    
    
    # top features
    tfs = pd.read_csv(sys.argv[2], index_col=0).iloc[:10, :]
    print(tfs)
    
    # others 
    halting_degree = int(sys.argv[3])
    out_dir = sys.argv[4]
    nodes_lookup = defaultdict(list)
    for e in parse_seq.parse_gfa(open(copangraph)):
        if e.type != parse_seq.GFATypes.S:
            continue
        nodes_lookup[e.nid].append(e)
    # loop through each important feature and extract graph
    extract_nodes = lambda s: re.findall('([0-9]+).*-> ([0-9]+)', s)[0]
    for f in tfs['feature']:
        
        # node 
        if f.startswith('n'):
            f = int(f[1:])
            start_node = nmap.loc[f, 0]
            print('start_node: ', start_node)
            subgraph_nodes = set()
            dfs(adjM, subgraph_nodes, start_node, halting_degree)
            bn = os.path.basename(sys.argv[2]).replace('.tf.csv', '')
            write_gfa_subgraph(adjM, nodes_lookup, edge_lookup, subgraph_nodes, os.path.join(OUT_DIR, f'{bn}_{start_node}.gfa'))
            write_contigs(sample_names, nodes_lookup, subgraph_nodes, os.path.join(OUT_DIR, f'{bn}_{start_node}.fasta'))
        
        # edge
        else:
            f = int(f[1:])
            start_edge = emap.loc[f, 0]
            print('start_edge: ', start_edge)
            u, v = extract_nodes(start_edge)
            start_edge = start_edge.replace(' -> ', '_').replace('(', '').replace(')', '')
            u_sg, v_sg = set(), set()
            dfs(adjM, u_sg, u, halting_degree)
            dfs(adjM, v_sg, v, halting_degree)
            bn = os.path.basename(sys.argv[2]).replace('.tf.csv', '')
            subgraph_nodes = u_sg if len(u_sg) > len(v_sg) else v_sg
            write_gfa_subgraph(adjM, nodes_lookup, edge_lookup, subgraph_nodes, os.path.join(OUT_DIR, f'{bn}_{start_edge}.gfa'))
            write_contigs(sample_names, nodes_lookup, subgraph_nodes, os.path.join(OUT_DIR, f'{bn}_{start_edge}.fasta'))