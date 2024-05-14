import os
import tqdm
import sys
import yaml
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from scipy import sparse as sp
from Utils import parse_seq as ps
from Utils.align.align.calign import aligner
from Utils.align.align.matrix import DNAFULL
import random
random.seed(42)
segments = None

def get_sequence_to_node_matrix(segments, nodes, min_overlap):
    # populate s2v matrix
    s2v_mat = sp.lil_matrix((len(segments), len(nodes)), dtype=np.uint8)
    for i, s in segments.items():
        if len(s.seq) >= min_overlap:
            s2v_mat[i, nodes[s.nid]] = 1
    s2v_mat = s2v_mat.tocsc()
    return s2v_mat

def reverse_complement(s):
    rcd = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return ''.join(
        rcd[e] for e in s
    )[::-1]

def palign(coords):
    global segments
    v, s1, s2 = coords

    # compute forward complement alignment
    s1_ori = segments[s1].rest[-1]
    s2_ori = segments[s2].rest[-1]
    if s1_ori == s2_ori:
        aln_obj = sorted(aligner(segments[s1].seq, segments[s2].seq, method='glocal', matrix=DNAFULL), key=lambda x: x.score, reverse=True)[0]
        edit_dist = 1 - (aln_obj.n_mismatches + aln_obj.n_gaps1 + aln_obj.n_gaps2) / len(aln_obj.seq1)
    else:
        # compute reverse complement alignment
        aln_obj = sorted(aligner(reverse_complement(segments[s1].seq), segments[s2].seq, method='glocal', matrix=DNAFULL), key=lambda x: x.score, reverse=True)[0]
        edit_dist = 1 - (aln_obj.n_mismatches + aln_obj.n_gaps1 + aln_obj.n_gaps2) / len(aln_obj.seq1)
    #return edit_dist, aln_obj, len(aln_obj.seq1)
    return v, edit_dist

def compute_comparison(s2v_mat, v):
    if (s2v_mat[:, v] > 0).sum() <= 1:
        return list()
    m = sp.tril(s2v_mat[:, v] @ s2v_mat[:, v].T, k=-1)
    if m.nnz == 0:
        return list()
    l_s1, l_s2, _ = sp.find(m)
    return list(zip(l_s1, l_s2))

def compute_within_node_similarity(segments, s2v_mat, sample_sz=10000, sample_proportion=None):
    
    # compute the alignments
    valid_nodes = list()
    print('Compute valid nodes...')
    for v in tqdm.tqdm(range(s2v_mat.shape[1])):
        if (s2v_mat[:, v] > 0).sum() <= 1:
            continue
        m = sp.tril(s2v_mat[:, v] @ s2v_mat[:, v].T, k=-1)
        valid_nodes.append(v)
    print('num_valid_nodes', len(valid_nodes))
        
    comparisons = list()
    print('Compute comparisons...')
    if sample_proportion:
        sample = random.sample(valid_nodes, int(len(valid_nodes) * sample_proportion))
    else:
        sample = random.sample(valid_nodes, sample_sz)
    for v in tqdm.tqdm(sample):
        m = sp.tril(s2v_mat[:, v] @ s2v_mat[:, v].T, k=-1)
        l_s1, l_s2, _ = sp.find(m)
        for s1, s2 in zip(l_s1, l_s2):
            comparisons.append((v, s1, s2))
    print('computing alignments between', sample_sz, 'nodes: ', len(comparisons))
    with Pool(processes=32) as pool:
        ret = pool.map(palign, comparisons)
        
    mean_aln_score_pv = np.zeros(len(sample))
    min_aln_score_pv = np.zeros(len(sample))
    ret_lookup = defaultdict(list)
    for node, ed in ret:
        ret_lookup[node].append(ed)
    for i, node_dists in enumerate(ret_lookup.values()): 
        min_aln_score_pv[i] = min(node_dists)
        mean_aln_score_pv[i] = np.mean(node_dists)
    return min_aln_score_pv, mean_aln_score_pv
    
        
def compute_between_node_similarity(segments, s2v_mat, min_overlap, sample_sz=10000):

    def get_valid_seq(seqs, mo):
        global segments 
        for i in seqs:
            if len(segments[i].seq) >= mo:
                return i
        return None

    num_nodes = s2v_mat.shape[1]
    # construct a set of valid sequence comparisons from different vertices. No two comparisons
    # can involve the same nodes.
    valid_seqs = list()
    node_pairs = set()
    print('selecting random sequences to align...', flush=True)
    while len(valid_seqs) != sample_sz:
        u, v = np.random.randint(0, num_nodes), np.random.randint(0, num_nodes)
        if u == v:
            continue

        if (u, v) in node_pairs or (v, u) in node_pairs:
            continue

        if s2v_mat[:, u].nnz == 0 or s2v_mat[:, v].nnz == 0:
            continue
        
        # get a sequence from each node longer than min_overlap, if one doesn't exist in each node, continue
        valid_u_seq = get_valid_seq(sp.find(s2v_mat[:, u])[0], min_overlap)
        valid_v_seq = get_valid_seq(sp.find(s2v_mat[:, v])[0], min_overlap)
        assert (valid_u_seq is not None and valid_v_seq is not None)
        
        # add seqs 
        node_pairs.add((u, v))
        node_pairs.add((v, u))
        valid_seqs.append((None, valid_u_seq, valid_v_seq))
        
    print('computing between node comparisons...', flush=True)
    with Pool(processes=32) as pool:
        ret = pool.map(palign, valid_seqs)

    scores = np.zeros(sample_sz)
    for i, (_, score) in enumerate(ret):
        scores[i] = score
    return scores


if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        print('Usage: <exe> <config>')
        sys.exit()
        
    # load config
    with open(sys.argv[1], 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        dn = config['dataset_name']
        gfa = config['gfa']
        out_dir = config['out_dir']
        min_overlap = int(config['min_overlap'])
    
    # load in gfa segments
    with open(gfa) as f:
        segments = [e for e in ps.parse_gfa(f) if e.type == ps.GFATypes.S]
    # construct node and segment lookups
    segments = {i:s for i, s in enumerate(segments)}
    nodes = {e.nid for e in segments.values()}
    nodes = {n:i for i, n in enumerate(sorted(list(nodes)))}
    
    # construct sequence to vertex matrix.
    # filter vertex include those with a sequence of len >= min_overlap
    s2n_mat = get_sequence_to_node_matrix(segments, nodes, min_overlap)
    
    # construct within node distances 
    within_node_min_similarity, within_node_mean_similarity = compute_within_node_similarity(segments, s2n_mat)
    
    # construct between sequence alignments
    between_node_similarity = compute_between_node_similarity(segments, s2n_mat, min_overlap)
    print('within - mean:', within_node_mean_similarity.mean())
    print('within - min:', within_node_min_similarity.mean())
    print('between- mean:', between_node_similarity.mean())
    
    # write to disk 
    np.save(os.path.join(out_dir, f'{dn}_within_node_min_similarity.npy'), within_node_min_similarity)
    np.save(os.path.join(out_dir, f'{dn}_within_node_mean_similarity.npy'), within_node_mean_similarity)
    np.save(os.path.join(out_dir, f'{dn}_between_node_similarity.npy'), between_node_similarity)