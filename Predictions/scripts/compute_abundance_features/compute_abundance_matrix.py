
import os
import sys
from tqdm import tqdm
import pandas as pd
import numpy as np
import parse_seq as ps

IDXSTATS = '/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/EXT/read_mappings'
SAMPLES = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/moms-pi/random.sample_names.csv'
GFA = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/moms-pi/momspi_sd005_edgesync2.gfa'
OUTFILE = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/moms-pi/momspi_sd005_edgesync2.abundance.csv'

def get_idx_df(samples):
    """Return a df of [sample name, contig name, contig len, depth]"""
    idx_dfs = list()
    library_size = dict()
    for e in samples:
        print('processing ', e)
        idxstat_file = f'{IDXSTATS}/{e}/mapping.coordsort.idxstat'
        df = pd.read_csv(idxstat_file, delimiter='\t', header=None)
        df.columns = ['contig_name', 'contig_len', 'depth', 'mate_mapped']
        df['sample_name'] = e
        library_size[e] = df.depth.sum()
        df.depth = (df.depth / df.contig_len) * 1000 # counts per kbp
        df.depth = df.depth / (library_size[e]/1e6) # control for library size
        df = df.loc[df.contig_name != '*',:]
        idx_dfs.append(df)
    idx_df = pd.concat(idx_dfs)
    return idx_df

def get_nodes(gfa_file):
    with open(gfa_file) as f:
        nodes = {r.nid for r in ps.parse_gfa(f) if r.type == ps.GFATypes.S}

    # make sure copangraph nodes are numbers 1:N
    for i in range(1, len(nodes)+1):
        assert str(i) in nodes
    return nodes

if __name__ == '__main__':
    with open(SAMPLES) as f:
        samples = [l.strip() for l in f]


    nodes = get_nodes(GFA)
    idx_df = get_idx_df(samples)
    idx_df.index = idx_df.apply(lambda x: f'{x.sample_name}:{x.contig_name}', axis=1)

    
    records = list()
    with open(GFA) as f:
        for r in tqdm(ps.parse_gfa(f)):
            if r.type != ps.GFATypes.S:
                continue
            records.append([r.nid, r.contig_name, r.sample_name, idx_df.loc[f'{r.sample_name}:{r.contig_name}', 'depth']])
    seq_depth_df = pd.DataFrame(records, columns=['nid', 'cn', 'sn', 'depth'])
    node_features = pd.DataFrame(np.full((len(samples), len(nodes)), 0.0), index=samples, columns =  sorted(list(nodes)))
    for (n,s), df in tqdm(seq_depth_df.groupby(['nid', 'sn'])):
        node_features.loc[s,n] = df.depth.sum()
    node_features.to_csv(OUTFILE)
