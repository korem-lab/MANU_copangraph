
import sys 
import os
import re
import pandas as pd
import numpy as np
import parse_seq 
from scipy import sparse as sp
from collections import defaultdict
from tqdm import tqdm

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
    for u, v, uori, vori in tqdm(edges):
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

def dfs(adjM, nodes, hairball_adjacent, dead_end, u, depth=5, hairball_degree=10, halt=False):
    u = str(u)
    nodes.add(u)
    if depth == 0:
        return
    

    adj_nodes = adjM[u]
    # hairball stoping condition
    if len(adj_nodes) >= hairball_degree:
        hairball_adjacent.add(u)
        return

    # tip removal 
    if len(adj_nodes) == 1:
        # remove the top and halt
        dead_end.add(u)
        nodes.remove(u)
        return

    # otherwise add nodes
    nodes.update(adj_nodes)
    for v in adj_nodes:
        dfs(adjM, nodes, hairball_adjacent, dead_end, v, depth-1, hairball_degree, halt)
        
def write_gfa_subgraph(adjM, nodes_lookup, edge_lookup, subgraph_nodes, gfa_name):
    with open(gfa_name, 'w') as f:
        for n in sorted(list(subgraph_nodes), key=lambda x: int(x)):
            samples = list(set([e.rest.split(':')[3] for e in nodes_lookup[n]]))
            depth = f'DP:f:{len(samples)}'
            for e in nodes_lookup[n]:
                e.depth = depth
                e.write(f)
        for u in sorted(list(subgraph_nodes), key=lambda x: int(x)):
            for v in subgraph_nodes:
                if v in adjM[u] or u in adjM[v]:
                    if u < v:
                        for _u, _v, uori, vori in edge_lookup[(u, v)]:
                            uori = '+' if uori==1 else '-'
                            vori = '+' if vori==1 else '-'
                            f.write(f'L\t{_u}\t{uori}\t{_v}\t{vori}\t0M\n')
                    else:
                        for _u, _v, uori, vori in edge_lookup[(v, u)]:
                            uori = '+' if uori==1 else '-'
                            vori = '+' if vori==1 else '-'
                            f.write(f'L\t{_u}\t{uori}\t{_v}\t{vori}\t0M\n')

def write_color_file(nodes_lookup, subgraph_nodes, outcome, filename):
    def red_blue_ratio(p,c):
        t = p+c
        frac = 255/t
        return int(np.floor(p * frac)), int(np.ceil(c * frac))

    def red_blue_grey(a, b, max_val):
        if a + b == 0:
            return '#d0d0d0'  # No color data, full white
    
        # Normalize red and blue contribution
        total = a + b
        red_ratio = a / total
        blue_ratio = b / total
    
        # Get base RGB from red-blue blend (no green in this case)
        r = red_ratio * 255
        g = 0
        b = blue_ratio * 255
        
        grey = 208
        # Whiteness factor: how far a+b is from max_val
        greyness = 1 - min(total / max_val, 1)
    
        # Blend towards grey
        r = r + (grey - r) * greyness
        g = g + (grey - g) * greyness
        b = b + (grey - b) * greyness
        return '#{0:02x}{1:02x}{2:02x}'.format(int(r), int(g), int(b))

    records = list()    
    for n in subgraph_nodes:
        samples = list(set([e.rest.split(':')[3] for e in nodes_lookup[n]]))
        num_persist = sum(outcome.loc[int(e), 'outcome'] for e in samples)
        num_clear = len(samples) - num_persist
        red_val, blue_val = red_blue_ratio(num_persist, num_clear)
        hex_color = '#' + '%.2X' % red_val + '00' + '%.2X' % blue_val
        records.append((n, hex_color, num_persist, num_clear))
    clr_df = pd.DataFrame(records, columns=['Node', 'ColourLocal', 'num_p', 'num_c'])
    max_samples = (clr_df.num_p + clr_df.num_c).max()
    clr_df['Colour'] = clr_df.apply(lambda x: red_blue_grey(x['num_p'], x['num_c'], max_samples), axis=1)
    clr_df = clr_df[['Node', 'Colour', 'ColourLocal', 'num_p', 'num_c']]
    clr_df.set_index('Node', inplace=True)
    clr_df.to_csv(filename)


def write_contigs(nodes_lookup, subgraph_nodes, fname, contig_dir):
    
    contigs = list()
    get_sample = lambda x: int(x.rest.split(':')[3])
    get_cn = lambda x: x.rest.split(':')[2]
    
    # get segments grouped by sample
    segments = defaultdict(list)
    for n in subgraph_nodes:
        for s in nodes_lookup[n]:
            segments[get_sample(s)].append(s)
            
    print('num_samples: ', len(segments)) 
    for sample_name, v in segments.items():
        cset = {get_cn(e) for e in v}
        with open(os.path.join(contig_dir, f'{sample_name}/final.contigs.fa')) as f:
            for fa in parse_seq.parse(f, parse_seq.Fasta):
                cn = fa.hdr.split(' ')[0]
                if cn in cset:
                    fa.hdr = f'{cn}_{sample_name}'
                    contigs.append(fa)
    print(len(contigs))
    with open(fname, 'w') as f:
        for c in contigs:
            c.write(f)

def write_bed_file(nodes_lookup, subgraph_nodes, fname):

    records = list()
    for n in subgraph_nodes:
        for s in nodes_lookup[n]:
            _, _, cn, sn, _, l, r, ori = s.rest.split(':')
            records.append(f'{sn}_{cn}\t{l}\t{r}\t{s.nid}')
    with open(fname, 'w') as f:
        f.write('\n'.join(records))

def write_hairball_adjacent(nodes, outfile):
    with open(outfile, 'w') as f:
        f.write('\n'.join(str(e) for e in nodes))
def write_deadend_nodes(nodes, outfile):
    with open(outfile, 'w') as f:
        f.write('\n'.join(str(e) for e in nodes))
            

if __name__ == '__main__':

    if len(sys.argv) != 7:
        print('Usage: <exe> <copangraph> <feature_list> <depth> <out_dir> <contig_dir> <outcome>')
        sys.exit()
    
    # load gfa 
    links = list()
    gfa, features, depth, outdir, contig_dir, outcome_fl = sys.argv[1:]
    outcome = pd.read_csv(outcome_fl, index_col=0)
    print(outcome.index)
    os.mkdir(outdir)

    depth = int(depth)

    graph_name = os.path.basename(gfa).replace('.gfa', '')
    adjM, edge_lookup = gfa2adjM(gfa)
    
    # others 
    nodes_lookup = defaultdict(list)
    with open(gfa) as _g:
        for e in parse_seq.parse_gfa(_g):
            if e.type == parse_seq.GFATypes.S:
                nodes_lookup[e.nid].append(e)

    # loop through each important feature and extract graph

    extract_nodes = lambda s: re.findall('([0-9]+).*-> ([0-9]+)', s)[0]


    features = pd.read_csv(features, header=None)[0].astype(str)
    for f in features:
        print(f)
        # node 
        if '->' not in f:
            print('start_node: ', f)
            subgraph_nodes = set()
            hairball_adjacent = set()
            dead_end= set()
            dfs(adjM, subgraph_nodes, hairball_adjacent, dead_end, f, depth)
            write_gfa_subgraph(adjM, nodes_lookup, edge_lookup, subgraph_nodes, os.path.join(outdir, f'{graph_name}_{f}.gfa'))
            write_contigs(nodes_lookup, subgraph_nodes, os.path.join(outdir, f'{graph_name}_{f}.fasta'), contig_dir)
            write_color_file(nodes_lookup, subgraph_nodes, outcome, os.path.join(outdir, f'{graph_name}_{f}.color.csv'))
            write_bed_file(nodes_lookup, subgraph_nodes, os.path.join(outdir, f'{graph_name}_{f}.subgraph.bed'))
            write_hairball_adjacent(hairball_adjacent, os.path.join(outdir, f'{graph_name}_{f}.hairball_nodes.txt'))
            write_deadend_nodes(dead_end, os.path.join(outdir, f'{graph_name}_{f}.deadend_nodes.txt'))
        
        # edge
        else:
            print('start_edge: ', start_edge)
            u, v = extract_nodes(start_edge)
            start_edge = start_edge.replace(' -> ', '_').replace('(', '').replace(')', '')
            u_sg, v_sg = set(), set()
            hairball_adjacent = set()
            dfs(adjM, u_sg, hairball_adjacent, u, halting_degree)
            dfs(adjM, v_sg, hairball_adjacent, v, halting_degree)
            bn = os.path.basename(sys.argv[2]).replace('.tf.csv', '')
            subgraph_nodes = u_sg | v_sg
            write_gfa_subgraph(adjM, nodes_lookup, edge_lookup, subgraph_nodes, os.path.join(outdir, f'{graph_name}_{f}.gfa'))
            write_contigs(nodes_lookup, subgraph_nodes, os.path.join(outdir, f'{graph_name}_{f}.fasta'))
