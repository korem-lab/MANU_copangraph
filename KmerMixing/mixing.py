import parse_seq as ps
import sys
import os
import glob
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import pandas as pd
import re
import random

K=21

def lex_low(k):
    rc = ''.join([{'A':'T', 'C':'G', 'G':'C', 'T':'A'}[e] for e in k][::-1])
    return k if k < rc else rc

def rev_comp(s):
    return ''.join([{'A':'T', 'C':'G', 'G':'C', 'T':'A'}[e] for e in s][::-1])

def draw_random_kmer(s, K=15):
    draw = random.randint(0, len(s)-K)
    return lex_low(s[draw: draw+K])


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

def get_kmer_sets(genome_list, K=15):
    kmer_sets = list()
    for i, g in tqdm(enumerate(genome_list)):
        kmer_sets.append(set())
        with open(g) as f:
            for fa in ps.parse(f, ps.Fasta):
                for j in range(len(fa.seq)-K+1):
                    kmer = lex_low(fa.seq[j: j+K])
                    kmer_sets[i].add(kmer)
    return kmer_sets

def parse_graph_to_node2seq(graph_file):
    if graph_file.endswith('.gfa'):
        is_copan = True
    elif graph_file.endswith('.fastg'):
        is_copan = False
    else:
        AssertionError
    
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
            n2seq[key].append(fa)
        else:
            # fastq files output by megahit_toolkit have 
            # each node duplicated in reverse complement
            # therefore, to avoid overcounting the nodes
            # remove one of the two (either forward or reverse
            # it doesn't matter for calculation)
            key = re.findall('(NODE_[0-9]+)_', fa.hdr).pop(0)
            if key not in n2seq:
                n2seq[key].append(fa)
    return n2seq

def seq_in_genome_mixing(graph_file, genome_files, tool, kmer_sets, candidate_its=1):
    genomes = list()
    for fl in genome_files:
        with open(fl) as f:
            for e in ps.parse(f, ps.Fasta):
                genomes.append(e)
    print('Num genomes read: ', len(genomes))
    if kmer_sets is None:
        print('Building kmer sets...', len(genomes))
        kmer_sets = get_kmer_sets(genome_files)
    print('parsing graph...')
    n2seq = parse_graph_to_node2seq(graph_file)
    data = list()
    for k, v in tqdm(n2seq.items()):
        max_bp = -1
        total_bp = 0
        unmapped_bp = 0
        genomes_in = set()
        for e in v:
            sequence_unmapped=True
            for i, g in enumerate(genomes):
                it = 0
                candidate = True
                while it < candidate_its:
                    if draw_random_kmer(e.seq) not in kmer_sets[i]:
                        candidate = False
                    if not candidate: 
                        break
                    it += 1
                if candidate:
                    mapped = (g.seq.find(e.seq) != -1 or g.seq.find(rev_comp(e.seq)) != -1)
                    if mapped:
                        genomes_in.add(i)
                        sequence_unmapped=False
            if sequence_unmapped:
                unmapped_bp += len(e.seq)
        for e in v:
            if len(e.seq) > max_bp:
                max_bp = len(e.seq)
            total_bp += len(e.seq)
        data.append([k, len(v), total_bp, max_bp, len(genomes_in) > 1,  len(genomes_in) == 1, unmapped_bp>0, unmapped_bp, len(genomes_in), tool])
    return pd.DataFrame(data, columns = [
            'nid', 'num_seqs', 'total_bp', 'max_bp', 'is_multi_genome', 'is_single_genome', 
            'has_unmapped', 'unmapped_bp', 'num_genomes_in', 'tool'
            ])


def kmer_to_nodes(graph_file):
    n2seq = parse_graph_to_node2seq(graph_file)
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

def n50(graph_file):

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
    data = sorted(data, key=lambda x: len(x.seq), reverse=True)
    total_50 = sum(len(e.seq) for e in data)/2
    cum_sum= 0 
    for e in data:
        cum_sum += len(e.seq)
        if cum_sum >= total_50:
            return len(e.seq)

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

