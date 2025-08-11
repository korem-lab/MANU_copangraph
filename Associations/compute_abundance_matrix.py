import os
import sys
import parse_seq as ps
import pandas as pd
import gzip
import tqdm
import numpy as np
import yaml

def get_num_nodes(gfa_file):
    with open(gfa_file) as f:
        nodes = {r.nid for r in ps.parse_gfa(f) if r.type == ps.GFATypes.S}

    # make sure copangraph nodes are numbers 1:N
    for i in range(1, len(nodes)+1):
        if str(i) not in nodes:
            print(i)
    return len(nodes)



def compute_abnd_mat(gfa_file, idx_stats_dir, map_samples, num_reads_series, short=500):
    """
     gfa_file: a copangraph made from samples in the graph_samples list
    """
   
    rbkms = dict()
    short_cntg = dict()
    for s in map_samples:
        idx_counts = pd.read_csv(os.path.join(idx_stats_dir, f'{s}_idxstats.csv'), delimiter='\t', header=None, index_col=0)[[2,3]].sum(axis=1)
        cntg_lns = pd.read_csv(os.path.join(idx_stats_dir, f'{s}_idxstats.csv'), delimiter='\t', header=None, index_col=0)[1]
        rbkm = idx_counts / ((cntg_lns/1e3) * (num_reads_series[s]/1e6))
        rbkms[s] = rbkm
        short_cntg[s] = cntg_lns < short
    
    # construct node features
    num_nodes = get_num_nodes(gfa_file)
    node_features = np.full((len(map_samples), num_nodes), 0.0)
    stoi = {s: i for i, s in enumerate(map_samples)}
    short = np.full((num_nodes,), False)
    count = 0
    with open(gfa_file) as f:
        for r in tqdm.tqdm(ps.parse_gfa(f)):
            if r.type != ps.GFATypes.S:
                continue
            # add to feature 
            s = r.sample_name
            c = r.contig_name
            rbkm = rbkms[s][c]
            node_features[stoi[s], int(r.nid)-1] += rbkm
            count += rbkm
        
            # record if short
            if short_cntg[s][c]:
                short[int(r.nid)-1] = True

    node_features = pd.DataFrame(node_features, index=map_samples, columns=range(1, num_nodes+1))
    # remove nodes deriving from short contigs
    return node_features.T, node_features.loc[:, ~short].T


if __name__ == '__main__':

    yml = sys.argv[1]
    assert(os.path.isfile(yml))
    with open(yml) as f:
        yml = yaml.safe_load(f)
    map_samples = pd.read_csv(yml['map_samples'], header=None, dtype=str)[0]
    idx_stats_dir = yml['idx_stats_dir']
    short_contig = yml['short_contig']
    num_reads = pd.read_csv(yml['read_counts'], index_col=0, header=None)
    num_reads = num_reads.loc[map_samples, 1]
    num_reads = num_reads.astype(float)
    gfa_file = yml['gfa_file']
    nf, len_filtered = compute_abnd_mat(gfa_file, idx_stats_dir, map_samples, num_reads, short_contig)
    nf.to_pickle(os.path.join(idx_stats_dir, os.path.basename(gfa_file).replace('.gfa', '') + '.abnd_mat.pkl'))
    len_filtered.to_pickle(os.path.join(idx_stats_dir, os.path.basename(gfa_file).replace('.gfa', '') + '.abnd_mat_filtered.pkl'))
    
