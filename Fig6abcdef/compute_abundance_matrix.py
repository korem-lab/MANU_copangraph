import os
import sys
import parse_seq_copangfa as ps
import pandas as pd
import gzip
import numpy as np
import yaml
import tqdm

def get_num_nodes(gfa_file):
    with open(gfa_file) as f:
        nodes = {r.nid for r in ps.parse_gfa(f) if r.type == ps.GFATypes.S}

    # make sure copangraph nodes are numbers 1:N
    for i in range(1, len(nodes)+1):
        if str(i) not in nodes:
            print(i)
    return len(nodes)



def compute_abnd_mat(gfa_file, idx_stats_dir, samples, num_reads, short=500):
    """
     gfa_file: a copangraph made from samples in the graph_samples list
    """
   
    rbkms = dict()
    short_cntg = dict()
    print('loading idxstats')
    for s in samples:
        idx_counts = pd.read_csv(os.path.join(idx_stats_dir, f'{s}.stats.csv'), delimiter='\t', header=None, index_col=0)[2]
        cntg_lns = pd.read_csv(os.path.join(idx_stats_dir, f'{s}.stats.csv'), delimiter='\t', header=None, index_col=0)[1]
        rbkm = idx_counts / ((cntg_lns/1e3) * (num_reads[s]/1e6))
        rbkms[s] = rbkm
        short_cntg[s] = cntg_lns < short
    
    # construct node features
    print('counting nodes')
    num_nodes = get_num_nodes(gfa_file)
    node_features = np.full((len(samples), num_nodes), 0.0, dtype=np.float32)
    print('node features', node_features.shape)
    print('making feature matrix')
    stoi = {s: i for i, s in enumerate(samples)}
    short = np.full((num_nodes,), False)
    with open(gfa_file) as f:
        for r in tqdm.tqdm(ps.parse_gfa(f)):
            if r.type != ps.GFATypes.S:
                continue
            # add to feature 
            s = r.sample_name
            c = r.contig_name
            rbkm = rbkms[s][c]
            node_features[stoi[s], int(r.nid)-1] = rbkm
            # record if short
            if short_cntg[s][c]:
                short[int(r.nid)-1] = True

    node_features = pd.DataFrame(node_features, index=samples, columns=range(1, num_nodes+1))
    # remove nodes deriving from short contigs
    return node_features.T, node_features.loc[:, ~short].T


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('exe: samples, num_reads, gfa_file, idx_stats_dir')
        sys.exit()
    samples, num_reads, gfa_file, idx_stats_dir = sys.argv[1:]
    samples = pd.read_csv(samples, header=None, dtype=str)[0]
    num_reads = pd.read_csv(num_reads, index_col=0, header=None)
    print(num_reads)
    num_reads = num_reads.loc[samples, 1]
    num_reads = num_reads.astype(float)
    abund, abund_filt = compute_abnd_mat(gfa_file, idx_stats_dir, samples, num_reads)
    abund.to_pickle(os.path.join(gfa_file.replace('.gfa', '') + '.abnd_mat.pkl'))
    
