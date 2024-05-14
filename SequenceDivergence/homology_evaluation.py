import os
import random

import tqdm
import pandas as pd
import pickle
import sys
import yaml
import numpy as np
from multiprocessing import Pool
from functools import partial
from scipy import sparse as sp
from scipy.stats import entropy
from align.calign import aligner
from align.matrix import DNAFULL

import src.parse_seq as ps
from evaluation.utils.manitou_io import monkey_patch_writer as mpw
from evaluation.utils.manitou_io import monkey_patch_reader as mpr

wopen = mpw(open)
ropen = mpr(open)

fastas = None
taxonomic_levels = ['Kingdom', 'Phylum', 'Order', 'Class', 'Family', 'Genera', 'Species', 'Strain']


def get_v2g(rg_fasta_list, g2c_mat, contig_map, vertex_map):
    """
    Constructs a vertex by genome matrix, where entries count the
    number of sequences in a vertex that come from a particular genome.

    A vertex by contig matrix is also produced as a by-product. And maps
    which map each named vertex or contig to the entry position in the
    numpy array.
    """

    # vertex x contig matrix: 1 if contig in vertex, 0 otherwise
    v2c_mat = sp.lil_matrix((len(vertex_map), len(contig_map)), dtype=np.uint32)
    for e in rg_fasta_list:
        v2c_mat[vertex_map[e.nid], contig_map[e.cn]] += 1

    # construct the vertex to genome map: V x G  = V x C @ C x G
    v2c_mat = sp.csc_matrix(v2c_mat)
    v2g_mat = v2c_mat @ g2c_mat.T
    return v2g_mat


def get_v2taxdstr(v2g_mat, g2tax):
    """
    Contstructs set of vertex by taxonomic distribution matrices that count the number
    of sequences (or species) that occurs in a node that belong to a particular taxa
    at a given taxonomic level.
    """
    # construct the taxonomy instance tables
    g2taxdstr_mats = list()
    g2taxdstr_maps = list()
    for taxa_lvl in taxonomic_levels:

        # for the given taxonomic level (e.g genera),
        # get the set of the genomes' assignments at that taxonomic level
        # (e.g the set of genera assignments of the genomes... {Escherichia, Salmonella, Lactobaccilus, ...}
        dstr_map = set(g2tax.loc[:, taxa_lvl])
        dstr_map = {e: i for i, e in enumerate(dstr_map)}

        # construct a matrix, G x (number of instances of a particular taxonomic level)
        # If, in the total set of genomes there are 5 genera, then the column number would be 5
        # Each entry in this matrix is 1 if the genome belongs to that taxonomic instance
        # and 0 otherwise
        g2taxinst_mat = sp.lil_matrix((g2tax.shape[0], len(dstr_map)), dtype=np.uint8)

        for g in g2tax.index:
            g2taxinst_mat[g, dstr_map[g2tax.loc[g, taxa_lvl]]] = 1

        # record the info in the dictionaries
        g2taxdstr_mats.append(sp.csc_matrix(g2taxinst_mat))
        g2taxdstr_maps.append(dstr_map)

    # construct the vertex to taxonomic distributions
    v2taxdstr_list = list()
    for taxa_lvl in range(len(taxonomic_levels)):
        g2taxdstr = g2taxdstr_mats[taxa_lvl]

        # V x (instances) = V x G @ G x (instances)
        v2taxdstr = v2g_mat @ g2taxdstr
        v2taxdstr_list.append(v2taxdstr)

    return v2taxdstr_list, g2taxdstr_mats, g2taxdstr_maps


def compute_entropy(v2taxdstr_list, diversity=False):
    entropies = list()
    for i, mat in enumerate(v2taxdstr_list):
        H = np.nan_to_num(entropy(mat.toarray(), axis=1))
        entropies.append(H if not diversity else np.exp(H))
    return entropies


def compute_genome_ANI(v2g, g2taxdstr_mats, ani_dists, avoid_self_comparison=True):
    k = -1 if avoid_self_comparison else 0
    # This function only makes sense with v2g as presence/absence
    v2g = v2g > 0
    comps = list()
    comps_avpv = list()
    for tax_lvl, g2taxdstr in enumerate(g2taxdstr_mats):
        print('Computing ANI at tax_level:', tax_lvl)

        # for each taxa at a given taxonomic level,
        # our goal is to construct a matrix G x G matrix that will be 1 for
        # genomes g,h if genomes g,h belong to the same taxa, 0 otherwise
        # for example, let taxa_lvl be genera: Lactobaccilus inners, and Lactobaccilus gasseri
        # would be 1, whilst L. inners and E. coli would be 0.
        #
        # To do this, we construt taxa_mat
        # - construct a rank 1 matrix from each taxa column
        # - multiply elementwise, this matrix by ani_dists
        # - tril (and 0 the diagonal if avoid_self_comparison == True)
        #
        # Once multiply this by ani_dists (element-wise) to get the matrix of distances
        # for this taxa
        num_comps = np.zeros(v2g.shape[0])
        total_dist = np.zeros(v2g.shape[0])
        comps = list()

        for taxa in range(g2taxdstr.shape[1]):
            # extract the vector denoting the each genomes occurrence in the current taxa
            g_taxa_occ_vec = g2taxdstr[:, taxa]

            # Multiply each row of V x G v2g matrix by this vector. This has the effect of
            # producing a matrix v2g_taxa where for each row, an entry is 1 if the corresponding
            # vertex contains that genome *and* that genome is within the taxa, 0 otherwise.
            # we convert to a numpy array to do a large matrix operation that follows
            v2g_taxa = v2g.multiply(g_taxa_occ_vec.T).toarray()

            # We construct a large V x G x G matrix, derived from v2g_taxa. For each row,
            # the G x G matrix is 1 if vertex v contains two genomes in the corresponding node *and*
            # both those genomes are from the taxa. To construct this matrix, we have to perform a broadcasting
            # multiplication, hence this set up. Essentially each G x G matrix is a rank-1 matrix
            # formed by multiplying each G row of the original V x G v2g_taxa by itself
            reshape_left = (v2g_taxa.shape[0], v2g_taxa.shape[1], 1)
            reshape_right = (v2g_taxa.shape[0], 1, v2g_taxa.shape[1])
            v2g_taxa = np.multiply(v2g_taxa.reshape(reshape_left), v2g_taxa.reshape(reshape_right))

            # We now multiply this by the ani distance matrix
            ani_dists = np.tril(ani_dists, k=k).astype('float32')
            ani_dists.reshape((1, ani_dists.shape[0], ani_dists.shape[1]))  # reshape to 1 x G x G
            v2g_taxa = np.multiply(v2g_taxa, ani_dists)

            # We can now calculate the per vertex comparisons and distance for this taxa
            num_comps += (v2g_taxa > 0).sum(axis=(1, 2))
            total_dist += v2g_taxa.sum(axis=(1, 2))
            v2g_taxa = v2g_taxa.flatten()
            comps.append(v2g_taxa[v2g_taxa != 0])

        # after calculating the distances, sum the total number of comparisons and distances to get the per vertex
        # and total average
        avg_pv = (np.divide(total_dist, num_comps, out=np.zeros_like(total_dist), where=num_comps != 0))
        comps_avpv.append(avg_pv)
        comps.append(np.hstack(comps))
    return comps_avpv, comps


def compute_genome_ANI_btwn(v2g, g2taxdstr_mats, ani_dists, avoid_self_comparison=True):
    k = -1 if avoid_self_comparison else 0
    # This function only makes sense with v2g as presence/absence
    v2g = v2g > 0
    comps = list()
    comps_avpv = list()

    # get all genomes present in a node
    reshape_left = (v2g.shape[0], v2g.shape[1], 1)
    reshape_right = (v2g.shape[0], 1, v2g.shape[1])
    v2g_all = np.multiply(v2g.toarray().reshape(reshape_left), v2g.toarray().reshape(reshape_right))
    v2g_all = v2g_all.astype(int)

    for tax_lvl, g2taxdstr in enumerate(g2taxdstr_mats):
        print('Computing ANI at tax_level:', tax_lvl)

        # for each taxa at a given taxonomic level,
        # our goal is to construct a matrix G x G matrix that will be 1 for
        # genomes g,h if genomes g,h belong to the same taxa, 0 otherwise
        # for example, let taxa_lvl be genera: Lactobaccilus inners, and Lactobaccilus gasseri
        # would be 1, whilst L. inners and E. coli would be 0.
        #
        # To do this, we construt taxa_mat
        # - construct a rank 1 matrix from each taxa column
        # - multiply elementwise, this matrix by ani_dists
        # - tril (and 0 the diagonal if avoid_self_comparison == True)
        #
        # Once multiply this by ani_dists (element-wise) to get the matrix of distances
        # for this taxa
        num_comps = np.zeros(v2g.shape[0])
        total_dist = np.zeros(v2g.shape[0])
        comps = list()
        v2g_within = np.zeros_like(v2g_all)

        for taxa in range(g2taxdstr.shape[1]):
            # extract the vector denoting the each genomes occurrence in the current taxa
            g_taxa_occ_vec = g2taxdstr[:, taxa]

            # Multiply each row of V x G v2g matrix by this vector. This has the effect of
            # producing a matrix v2g_taxa where for each row, an entry is 1 if the corresponding
            # vertex contains that genome *and* that genome is within the taxa, 0 otherwise.
            # we convert to a numpy array to do a large matrix operation that follows
            v2g_taxa = v2g.multiply(g_taxa_occ_vec.T).toarray()

            # We construct a large V x G x G matrix, derived from v2g_taxa. For each row,
            # the G x G matrix is 1 if vertex v contains two genomes in the corresponding node *and*
            # both those genomes are from the taxa. To construct this matrix, we have to perform a broadcasting
            # multiplication, hence this set up. Essentially each G x G matrix is a rank-1 matrix
            # formed by multiplying each G row of the original V x G v2g_taxa by itself
            reshape_left = (v2g_taxa.shape[0], v2g_taxa.shape[1], 1)
            reshape_right = (v2g_taxa.shape[0], 1, v2g_taxa.shape[1])
            v2g_taxa = np.multiply(v2g_taxa.reshape(reshape_left), v2g_taxa.reshape(reshape_right))

            # add to within
            v2g_within = v2g_within + v2g_taxa

        # take within taxa comparisons away from all comparisons to get between taxa comparisons
        v2g_btwn = v2g_all - v2g_within

        # We now multiply this by the ani distance matrix
        ani_dists = np.tril(ani_dists, k=k).astype('float32')
        ani_dists.reshape((1, ani_dists.shape[0], ani_dists.shape[1]))  # reshape to 1 x G x G
        v2g_btwn = np.multiply(v2g_btwn, ani_dists)

        # We can now calculate the per vertex comparisons and distance for this taxa
        num_comps += (v2g_btwn > 0).sum(axis=(1, 2))
        total_dist += v2g_btwn.sum(axis=(1, 2))
        v2g_btwn = v2g_btwn.flatten()
        comps.append(v2g_btwn[v2g_btwn!= 0])

        # after calculating the distances, sum the total number of comparisons and distances to get the per vertex
        # and total average
        avg_pv = (np.divide(total_dist, num_comps, out=np.zeros_like(total_dist), where=num_comps != 0))
        print(tax_lvl, avg_pv[avg_pv != 0].mean())
        comps_avpv.append(avg_pv)
        comps.append(np.hstack(comps))
    return comps_avpv, comps


def get_sequence_matrices(fas, g2c_mat, contig_map, vertex_map):
    # S x C
    s2c_mat = sp.lil_matrix((len(fas), len(contig_map)), dtype=np.uint8)
    # S x V
    s2v_mat = sp.lil_matrix((len(fas), len(vertex_map)), dtype=np.uint8)
    for i, fa in enumerate(fas):
        s2c_mat[i, contig_map[fa.cn]] = 1
        s2v_mat[i, vertex_map[fa.nid]] = 1
    s2c_mat = s2c_mat.tocsc()
    s2v_mat = s2v_mat.tocsc()

    # S x G
    s2g_mat = s2c_mat @ g2c_mat.T
    return s2c_mat, s2g_mat, s2v_mat


def reverse_complement(s):
    rcd = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return ''.join(
        rcd[e] for e in s
    )[::-1]


def palign(coords):
    global fastas
    s1, s2 = coords

    # compute forward complement alignment
    if fastas[s1].ori == fastas[s2].ori:
        aln_obj = sorted(aligner(fastas[s1].seq, fastas[s2].seq, method='glocal', matrix=DNAFULL), key=lambda x: x.score, reverse=True)[0]
        edit_dist = 1 - (aln_obj.n_mismatches + aln_obj.n_gaps1 + aln_obj.n_gaps2) / len(aln_obj.seq1)
    else:
        # compute reverse complement alignment
        aln_obj = sorted(aligner(reverse_complement(fastas[s1].seq), fastas[s2].seq, method='glocal', matrix=DNAFULL), key=lambda x: x.score, reverse=True)[0]
        edit_dist = 1 - (aln_obj.n_mismatches + aln_obj.n_gaps1 + aln_obj.n_gaps2) / len(aln_obj.seq1)
    return edit_dist, aln_obj, len(aln_obj.seq1)


def construct_alignment_matrix(s2v_mat, fas, out, dt, sample_id):

    s2v_mat = (s2v_mat != 0)
    aln_scores = sp.lil_matrix((s2v_mat.shape[0], s2v_mat.shape[0]))
    aln_lens = sp.lil_matrix((s2v_mat.shape[0], s2v_mat.shape[0]))

    # if alignment objects already exist, dont bother recomputing
    if os.path.exists(os.path.join(out, f'{sample_id}_alignment_objs_{dt}.pkl')):
        print('alignment objs already exist, loading...', flush=True)
        with wopen(os.path.join(out, f'{sample_id}_alignment_objs_{dt}.pkl'), 'rb') as f:
            alignment_objs = pickle.load(f)
        for s1, s2, aln_obj, aln_len in alignment_objs:
            edit_dist = 1 - (aln_obj.n_mismatches + aln_obj.n_gaps1 + aln_obj.n_gaps2) / len(aln_obj.seq1)
            aln_scores[s1, s2] = edit_dist
            aln_scores[s2, s1] = edit_dist
            aln_lens[s1, s2] = len(aln_obj.seq1)
            aln_lens[s2, s1] = len(aln_obj.seq1)
    else:
        print('no alignment objs found, constructing from scratch', flush=True)
        alignment_objs = list()
        with Pool(processes=16) as pool:

            for v in tqdm.tqdm(range(s2v_mat.shape[1])):
                # Get the pair-wise combinations of all sequence in a node. Compute an alignment
                # score for each
                if (s2v_mat[:, v] > 0).sum() <= 1:
                    continue
                m = sp.tril(s2v_mat[:, v] @ s2v_mat[:, v].T, k=-1)
                if m.nnz == 0:
                    continue
                l_s1, l_s2, _ = sp.find(m)

                ret = pool.map(palign, zip(l_s1, l_s2))

                for s1, s2, (edit_dist, aln_obj, aln_len) in zip(l_s1, l_s2, ret):
                    aln_scores[s1, s2] = edit_dist
                    aln_scores[s2, s1] = edit_dist
                    aln_lens[s1, s2] = aln_len
                    aln_lens[s2, s1] = aln_len
                    alignment_objs.append((s1, s2, aln_obj, aln_len))

        with wopen(os.path.join(out, f'{sample_id}_alignment_objs_{dt}.pkl'), 'wb') as f:
            pickle.dump(alignment_objs, f)

    return sp.tril(sp.csc_matrix(aln_scores), k=-1), sp.tril(sp.csc_matrix(aln_lens), k=-1)


def compute_alignment_scores(s2g_mat, s2v_mat, g2taxdstr_mats, alignment_scores, alignment_lengths=None):
    k = -1
    scores_total = list()
    scores_avpv = list()
    scores_minpv = list()
    alignment_scores = alignment_scores.tocsr()
    for tax_lvl, g2taxdstr in enumerate(g2taxdstr_mats):
        print('Computing alignment scores at tax_level:', tax_lvl, flush=True)
        if tax_lvl > 0:
            break

        # S x Taxa = S x G @ G x Taxa
        # For each sequence, get its taxonomic assignment
        s2taxdstr = (s2g_mat @ g2taxdstr != 0)

        scores_for_taxlvl = list()
        scores_pv = np.zeros(s2v_mat.shape[1])
        min_score_pv = np.full(s2v_mat.shape[1], np.inf)
        comps_pv = np.zeros(s2v_mat.shape[1])
        print("S2V SHAPE: ", s2v_mat.shape)
        for taxa in range(g2taxdstr.shape[1]):

            # Multiply each column of the sequence to vertex matrix by a
            # column  of the sequence to taxa matrix. This constructs a matrix that
            # is 1 if a sequence is in a vertex and is part of the taxonomic group,
            # and is 0 otherwise.
            s2v_taxa = s2v_mat.multiply(s2taxdstr[:, taxa]).tocsr()

            # For each vertex (column) retrieve the alignment scores
            # for each sequence within the vertex.
            for v in tqdm.tqdm(range(s2v_taxa.shape[1])):
                # retrieve the alignment scores

                # exit loop early
                if s2v_taxa[:, v].nnz <= 1:
                    continue

                occ_mat = s2v_taxa[:, v] @ s2v_taxa[:, v].T
                occ_mat = sp.tril(occ_mat, k=k)

                # exit loop early
                if occ_mat.nnz == 0:
                    continue

                score_mat = occ_mat.multiply(alignment_scores)
                if alignment_lengths is None:
                    score_mat_r, score_mat_c, dense_scores = sp.find(score_mat)
                    scores_pv[v] += score_mat.sum()
                    min_score_pv[v] = min(min_score_pv[v], dense_scores.min())
                    comps_pv[v] += score_mat.nnz
                else:
                    aln_len_mat = occ_mat.multiply(alignment_lengths)
                    score_mat = score_mat.multiply(aln_len_mat)
                    score_mat_r, score_mat_c, dense_scores = sp.find(score_mat)
                    scores_pv[v] += score_mat.sum()
                    comps_pv[v] += aln_len_mat.sum()

                    # get the index of the minimum score in the dense score array
                    min_score = float(dense_scores.min())
                    _, min_index, _ = sp.find(dense_scores == min_score)
                    min_len = aln_len_mat[score_mat_r[min_index[0]], score_mat_c[min_index[0]]]
                    assert min_score == score_mat[score_mat_r[min_index[0]], score_mat_c[min_index[0]]]
                    min_score_pv[v] = min(min_score_pv[v], min_score/min_len)
                scores_for_taxlvl.append(dense_scores)

        scores_total.append(np.hstack(scores_for_taxlvl))
        # compute the per-vertex average alignment score, and record this list in scores_avpv
        scores_minpv.append(min_score_pv)
        scores_avpv.append(np.divide(scores_pv, comps_pv, out=np.zeros_like(scores_pv), where=comps_pv != 0))

    return scores_total, scores_avpv, scores_minpv


def between_node_distances(s2v_mat, samples=1000, min_overlap=0):

    def get_valid_seq(seqs, mo):
        global fastas
        for i in seqs:
            if len(fastas[i].seq) >= mo:
                return i
        return None

    num_seqs, num_nodes = s2v_mat.shape
    # construct a set of valid sequence comparisons from different vertices. No two comparisons
    # can involve the same nodes.
    valid_seqs = list()
    node_pairs = set()
    print('Constructing randoms sequence pairs', flush=True)
    while len(valid_seqs) != samples:
        u, v = np.random.randint(0, num_nodes), np.random.randint(0, num_nodes)
        if u == v:
            continue

        if (u, v) in node_pairs or (v, u) in node_pairs:
            continue

        if s2v_mat[:, u].nnz == 0 or s2v_mat[:, v].nnz == 0:
            continue

        # for each node, make sure they have at least one sequence > min_overlap
        # since copan can only collapse sequences > min_overlap, only compare such sequences

        # get a sequence from each node longer than min_overlap, if one doesn't exist in each node, continue
        valid_u_seq = get_valid_seq(sp.find(s2v_mat[:, u])[0], min_overlap)
        valid_v_seq = get_valid_seq(sp.find(s2v_mat[:, v])[0], min_overlap)
        if (valid_u_seq is None) or (valid_v_seq is None):
            continue

        # if reached here, take the sequneces and record the node comparison
        node_pairs.add((u, v))
        node_pairs.add((v, u))
        valid_seqs.append((valid_u_seq, valid_v_seq))
    print('aligning selected sequence pairs', flush=True)
    with Pool(processes=16) as pool:
        ret = pool.map(palign, valid_seqs)

    scores = np.zeros(samples)
    for i, ((us, vs), (score, aln_obj, aln_len)) in enumerate(zip(valid_seqs, ret)):
        scores[i] = score
    return scores

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("usage: <exe> <config>")
        sys.exit(-1)

    with ropen(sys.argv[1], 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        dt = config['distance_threshold']
        nodes_fa = config['nodes_fasta']
        gc_occ = config['genome_contig_occurrence_matrix']
        gg_dist = config['genome_genome_distance_matrix']
        gtax = config['genome_taxonomy_matrix']
        out = config['out_directory']
        sample_id = config['sample_id']
        depth = config['depth']
        asm = config['asm']

    # read in g2c
    g2c_mat = pd.read_csv(gc_occ, index_col=0)
    genome_map = {g: i for i, g in enumerate(g2c_mat.index)}
    g2c_mat.index = range(g2c_mat.shape[0])

    # read in fasta elements
    fastas = list()
    with ropen(nodes_fa, 'r') as f:
        for e in ps.parse(f, ps.Fasta):
            fa = ps.CopanNodeSeq(*e.hdr.split(':'), seq=e.seq)
            if len(fa.seq) == 0:
                continue
            fastas.append(fa)


    # get a mapping from contigs to numpy index
    # note that, some contigs don't have a genome mapping measurement, so we must filter out the fasts
    # to only include those for which a measurement was obtained.
    contig_map = set(g2c_mat.columns)
    fastas = [e for e in fastas if e.cn in contig_map]
    g2c_mat = g2c_mat.loc[:, sorted(list(contig_map))]  # sort the contigs
    contig_map = {e: i for i, e in enumerate(sorted(list(contig_map)))}
    g2c_mat.columns = list(range(g2c_mat.shape[1]))
    g2c_mat = sp.csc_matrix(g2c_mat)
    print('num sequences after contig map filtering: ', len(fastas))

    # get a mapping from nodes to numpy index
    vertex_map = {e.nid for e in fastas}
    vertex_map = {e: i for i, e in enumerate(sorted(list(vertex_map)))}
    print('num sequences: ', len(fastas), ' num vertex ',len(vertex_map))
    
    # get genome to taxonomy table
    g2tax = pd.read_csv(gtax, index_col=0)
    g2tax.index = g2tax.file_name
    current_order = [genome_map[g] for g in g2tax.index]
    g2tax.index = current_order
    g2tax = g2tax.loc[list(range(g2tax.shape[0])), :]  # ensure g2tax is in the same order as the genome map
    taxa_lvl_map = {e: i for i, e in enumerate(g2tax.columns)}

    # parse ANI dists between genomes
    ani_dists = pd.read_csv(gg_dist, index_col=0, delimiter='\t')
    genomes = [e.replace('.fa', '') for e in genome_map.keys()]
    ani_dists = ani_dists.loc[genomes, genomes]  # filter to current genome selection
    ani_dists.index = [f'{e}.fa' for e in ani_dists.index]
    ani_dists.columns = [f'{e}.fa' for e in ani_dists.columns]
    current_order = [genome_map[g] for g in ani_dists.index]
    ani_dists.index = current_order
    ani_dists = ani_dists.loc[list(range(ani_dists.shape[0])), :]
    current_order = [genome_map[g] for g in ani_dists.columns]
    ani_dists.columns = current_order
    ani_dists = ani_dists.loc[:, list(range(ani_dists.shape[1]))]

    # compute matrices
    v2g_mat = get_v2g(fastas, g2c_mat, contig_map, vertex_map)
    v2taxdstr_list, g2taxdstr_mats, g2taxdstr_maps = get_v2taxdstr(v2g_mat, g2tax)

    # compute alignment
    print('Computing Alignments', flush=True)
    s2c_mat, s2g_mat, s2v_mat = get_sequence_matrices(fastas, g2c_mat, contig_map, vertex_map)
    
    alignment_scores, alignment_lengths = construct_alignment_matrix(s2v_mat, fastas, out, dt, sample_id)
    aln_scores, aln_scores_avpv, aln_scores_minpv = compute_alignment_scores(s2g_mat, s2v_mat, g2taxdstr_mats, alignment_scores, alignment_lengths=alignment_lengths)

    # compute between-node distances
    print('between node distances', flush=True)
    btwn_node_dists = between_node_distances(s2v_mat, samples=5000)

    with wopen(os.path.join(out, f'{sample_id}_btwn_seq_{dt}.pkl'), 'wb') as f:
        pickle.dump(btwn_node_dists, f)
    with wopen(os.path.join(out, f'{sample_id}_aln_{dt}.pkl'), 'wb') as f:
        pickle.dump(aln_scores, f)

    aln_scores_avpv = pd.DataFrame(aln_scores_avpv, index=taxonomic_levels)
    aln_scores_avpv['depth'] = depth
    aln_scores_avpv['sample_id'] = sample_id
    aln_scores_avpv['asm'] = asm
    aln_scores_avpv['dt'] = dt
    aln_scores_avpv.to_pickle = mpw(aln_scores_avpv.to_pickle)
    aln_scores_avpv.to_pickle(os.path.join(out, f'{sample_id}_aln_avpv_{dt}.pkl'))

    aln_scores_minpv = pd.DataFrame(aln_scores_minpv, index=taxonomic_levels)
    aln_scores_minpv['depth'] = depth
    aln_scores_minpv['sample_id'] = sample_id
    aln_scores_minpv['asm'] = asm
    aln_scores_minpv['dt'] = dt
    aln_scores_minpv.to_pickle = mpw(aln_scores_minpv.to_pickle)
    aln_scores_minpv.to_pickle(os.path.join(out, f'{sample_id}_aln_minpv_{dt}.pkl'))

    print('Computing Entropy', flush=True)
    entropies = compute_entropy(v2taxdstr_list)
    entropies = pd.DataFrame(entropies, index=taxonomic_levels)
    entropies['depth'] = depth
    entropies['sample_id'] = sample_id
    entropies['asm'] = asm
    entropies['dt'] = dt
    entropies.to_pickle = mpw(entropies.to_pickle)
    entropies.to_pickle(os.path.join(out, f'{sample_id}_entropies_{dt}.pkl'))

    #print('Computing ANI', flush=True)
    #ani_avpv, ani = compute_genome_ANI(v2g_mat, g2taxdstr_mats, ani_dists)
    #ani_avpv = pd.DataFrame(ani_avpv, index=taxonomic_levels)
    #ani_avpv['depth'] = depth
    #ani_avpv['sample_id'] = sample_id
    #ani_avpv['asm'] = asm
    #ani_avpv['dt'] = dt
    #ani_avpv.to_pickle = mpw(ani_avpv.to_pickle)
    #ani_avpv.to_pickle(os.path.join(out, f'{sample_id}_ani_avpv_{dt}.pkl'))
    #with open(os.path.join(out, f'{sample_id}_ani_{dt}.pkl'), 'wb') as f:
    #    pickle.dump(ani, f)

    #print('Computing btwn ANI', flush=True)
    #ani_avpv, ani = compute_genome_ANI_btwn(v2g_mat, g2taxdstr_mats, ani_dists)
    #ani_avpv = pd.DataFrame(ani_avpv, index=taxonomic_levels)
    #ani_avpv['depth'] = depth
    #ani_avpv['sample_id'] = sample_id
    #ani_avpv['asm'] = asm
    #ani_avpv['dt'] = dt
    #ani_avpv.to_pickle = mpw(ani_avpv.to_pickle)
    #ani_avpv.to_pickle(os.path.join(out, f'{sample_id}_btwn_ani_avpv_{dt}.pkl'))
    #with open(os.path.join(out, f'{sample_id}_btwn_ani_{dt}.pkl'), 'wb') as f:
    #    pickle.dump(ani, f)
