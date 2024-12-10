import parse_seq as ps
import sys
import os
import glob
import numpy as np
from collections import defaultdict
import pandas as pd
import re

K=21

def lex_low(k):
    rc = ''.join([{'A':'T', 'C':'G', 'G':'C', 'T':'A'}[e] for e in k][::-1])
    return k if k < rc else rc

def kmer_to_sample(sample_list, mode):

    k2s = defaultdict(int)
    for sample in sample_list:
        occurred_in_sample = set()
        with open(sample) as f:
            for fa in ps.parse(f, ps.Fasta):
                print(fa)
                for i in range(len(fa.seq)-K+1):
                    kmer = lex_low(fa.seq[i: i+K])
                    if mode == 'freq':
                        k2s[kmer] += 1
                    elif mode == 'occ' and (kmer not in occurred_in_sample):
                        k2s[kmer] += 1
                        occurred_in_sample.add(kmer)

    return k2s


def kmer_to_nodes(graph_file):

    if graph_file.endswith('.gfa'):
        is_copan = True
    else:
        is_copan = False
    
    with open(graph_file) as f:
        if is_copan:
            data = list(e for e in ps.parse_gfa(f) if e.type == ps.GFATypes.S)
        else:
            data = list(ps.parse(f, ps.Fasta))

    # organize sequences by node
    n2seq = defaultdict(list)
    for fa in data:
        if is_copan:
            key = fa.nid # node name
            n2seq[key].append(fa.seq)
        else:
            # fastq files output by megahit_toolkit have 
            # each node duplicated in reverse complement
            # therefore, to avoid overcounting the nodes
            # remove one of the two (either forward or reverse
            # it doesn't matter for calculation)
            key = re.findall('(NODE_[0-9]+)_', fa.hdr).pop(0)
            if key not in n2seq:
                n2seq[key].append(fa.seq)

    k2n = defaultdict(int) 
    for _,v in n2seq.items():
        in_node = set()
        for seq in v:
            for i in range(len(seq)-K+1):
                kmer = lex_low(seq[i:i+K])
                if kmer not in in_node:
                    k2n[kmer] += 1
                    in_node.add(kmer)
    return k2n 
    


def calc_mixing(num_samples, num_nodes):
    if num_samples == 1 and num_nodes == 1:
        return 1
    if num_samples == 1: # and num_nodes > 1
        return 1 - num_nodes
    return (float(num_samples) - num_nodes) / (float(num_samples) - 1)

def num_nodes(graph_file):

    if graph_file.endswith('.gfa'):
        is_copan = True
    else:
        is_copan = False
    
    with open(graph_file) as f:
        if is_copan:
            data = list(e for e in ps.parse_gfa(f) if e.type == ps.GFATypes.S)
        else:
            data = list(ps.parse(f, ps.Fasta))

    # organize sequences by node
    nodes = set()
    for fa in data:
        if is_copan:
            key = fa.nid
            nodes.add(key)
        else:
            # fastq files output by megahit_toolkit have 
            # each node duplicated in reverse complement
            # therefore, to avoid overcounting the nodes
            # remove one of the two (either forward or reverse
            # it doesn't matter for calculation)
            key = re.findall('(NODE_[0-9]+)_', fa.hdr).pop(0)
            nodes.add(key)

    return len(nodes)

def mixing(graph_file, k2s):


    # node lookup
    print('node calc...')
    k2n = kmer_to_nodes(graph_file)

    print('mixing calc...')
    missing = 0
    mixing_vec = list()
    for k in k2s.keys():
        if k2n[k] == 0:
            missing += 1
        else:
            m = calc_mixing(k2s[k], k2n[k])
            mixing_vec.append(m)

    mixing_vec = pd.Series(mixing_vec)
    return mixing_vec, float(missing)/len(k2s), len(k2s), len(k2n)

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print('Usage: <exe> <sample_dir> <seq_graph> <mode> <out_vec_name>')
        sys.exit()
    
    _, sample_dir, seqgraph_dir, graph, mode, outname = sys.argv
    assert(mode in ['occ', 'freq'])

    samples = glob.glob(os.path.join(sample_dir, '*.fasta'))
    mixing = pd.DataFrame(mixing(samples, graph, mode))
    mixing.to_pickle(outname)

