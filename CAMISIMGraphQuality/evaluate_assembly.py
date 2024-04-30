import numpy as np
import pandas as pd
from intervaltree import IntervalTree, Interval
from collections import defaultdict, namedtuple
from itertools import product, combinations
from sklearn.cluster import AgglomerativeClustering
from scipy import sparse as sp
import tqdm
from multiprocessing.dummy import Pool as TPool
from multiprocessing import Pool as PPool
import functools

Alignment = namedtuple('Alignment', ['contig_name', 'start', 'end', 'genome', 'contig_len', 'is_primary'])


def get_top_alignments(alignments, mode='is_primary'):
    assert mode in [
        'is_primary', 'mapq', 'AS'
    ]

    top_alignments = set()
    if mode == 'is_primary':
        for i, a in enumerate(alignments):
            if a.is_primary:
                top_alignments.add(i)
        return top_alignments

    if mode == 'mapq':
        score = lambda x: int(x.mapq)
    elif mode == 'AS':
        score = lambda x: float(x.AS)

    # break up by contig
    by_contig = defaultdict(list)
    for i, a in enumerate(alignments):
        by_contig[a.contig_name].append((i, a))

    # find the best alignment per contig and save it's index in alignments
    for k, v in by_contig.items():
        max_idx = -1
        max_score = -np.inf
        for i, a in v:
            if a.score > max_score:
                max_idx = i
                max_score = score(a)
        top_alignments.add(max_idx)
    return top_alignments


def get_coverage_metrics(asm_name, alignments, genomes, genome_lens, asm, g_to_aln_index, unassigned_contigs, contig_map, minimap_kmer_size=15):
    # get record objects
    tp_records = {
        g: {'genome': g, 'assembler': asm_name, 'metric': 'cov_tp', 'value': 0} for g in genomes
    }

    # we can have a per-genome fp, as sequences that "overcover" a genome, as well as non-genome-specific fp, which
    # fail to map anywhere
    fp_records = {
        g: {'genome': g, 'assembler': asm_name, 'metric': 'cov_fp', 'value': 0} for g in genomes
    }
    fp_records['-'] = {'genome': '-', 'assembler': asm_name, 'metric': 'cov_fp', 'value': 0}

    fn_records = {
        g: {'genome': g, 'assembler': asm_name, 'metric': 'cov_fn', 'value': 0} for g in genomes
    }

    #for i, g in tqdm.tqdm(enumerate(genomes)):
    #    aln_on_g = sorted(list(g_to_aln_index[g]))
    #    if len(aln_on_g) == 0:
    #        continue
    #    C = sp.lil_matrix((genome_lens[i], len(aln_on_g)), dtype=bool)
    #    Index = sp.lil_matrix((len(aln_on_g), len(contig_map)), dtype=bool)
    #    #cov_vector = np.zeros(genome_lens[i], dtype=np.uint32)
    #    #for i, aln_idx in enumerate(aln_on_g):
    #    #    a = alignments[aln_idx]
    #    #    cov_vector[a.start:a.end] += 1

    #    for j, aln_idx in enumerate(aln_on_g):
    #        a = alignments[aln_idx]
    #        C[a.start:a.end, j] = 1
    #        Index[j, contig_map[a.contig_name]] = 1
    #    C = sp.csr_matrix(C)
    #    Index = sp.csr_matrix(Index)
    #    Cov = C @ Index
    #    Cov = (Cov > 0)
    #    cov_vector = Cov.sum(axis=1)
    #    tp_records[g]['value'] = np.sum(cov_vector > 0)
    #    fn_records[g]['value'] = np.sum(cov_vector == 0)
    #    fp_records[g]['value'] = np.sum(cov_vector > 1)
    
    
    fun = functools.partial(get_genome_cov, g_to_aln_index, genome_lens, contig_map, alignments) 
    with PPool(processes=16) as pool:
        ret = pool.map(fun, enumerate(genomes))
        for dat in ret:
            if dat is None:
                continue
            g, tp, fn, fp = dat
            tp_records[g]['value'] = tp
            fn_records[g]['value'] = fn
            fp_records[g]['value'] = fp
        

    # contigs that fail to align to any genome
    num_unassigned_bp = 0
    for idx in unassigned_contigs:
        n, seq = asm.contigs[idx]
        if len(seq) < minimap_kmer_size:
            continue
        num_unassigned_bp += len(seq)
    fp_records['-']['value'] = num_unassigned_bp

    return list(tp_records.values()) + list(fn_records.values()) + list(fp_records.values())
        
def get_genome_cov(g_to_aln_index, genome_lens, contig_map, alignments, dat):
    i, g = dat
    aln_on_g = sorted(list(g_to_aln_index[g]))
    if len(aln_on_g) == 0:
        return None
    C = sp.lil_matrix((genome_lens[i], len(aln_on_g)), dtype=bool)
    Index = sp.lil_matrix((len(aln_on_g), len(contig_map)), dtype=bool)
    #cov_vector = np.zeros(genome_lens[i], dtype=np.uint32)
    #for i, aln_idx in enumerate(aln_on_g):
    #    a = alignments[aln_idx]
    #    cov_vector[a.start:a.end] += 1

    for j, aln_idx in enumerate(aln_on_g):
        a = alignments[aln_idx]
        C[a.start:a.end, j] = 1
        Index[j, contig_map[a.contig_name]] = 1
    C = sp.csr_matrix(C)
    Index = sp.csr_matrix(Index)
    Cov = C @ Index
    Cov = (Cov > 0)
    cov_vector = Cov.sum(axis=1)
    return (g, np.sum(cov_vector > 0), np.sum(cov_vector == 0), np.sum(cov_vector > 1))


def get_connectivity_metrics(
        asm_name,
        alignments,
        genomes,
        g_to_aln_idx,
        node_to_aln_idx,
        breakpoint_dict,
        asm,
        assembly_artifact_gap,
        window_sz,
        unaligned
):
    # get record objects
    tp_records = {
        g: {'genome': g, 'assembler': asm_name, 'metric': 'cnx_tp', 'value': 0} for g in genomes
    }
    fp_record = {'genome': '-', 'assembler': asm_name, 'metric': 'cnx_fp', 'value': 0}
    fn_records = {
        g: {'genome': g, 'assembler': asm_name, 'metric': 'cnx_fn', 'value': 0} for g in genomes
    }

    # get coordinates of all edges
    r, c, _ = sp.find(sp.tril(asm.adjM.sparse.to_coo()))

    # iterate edges and count FP
    print('Computing FP...')
    unaligned = {asm.contigs[i][0] for i in unaligned}
    for u, v in tqdm.tqdm(zip(asm.adjM.index[r], asm.adjM.index[c])):

        # If the edge involves an unaligned contig, we cannot measure its status
        if u in unaligned or v in unaligned:
           continue

        valid_edge = False
        print('\n\n\n',u, v)
        aln_found = False
        for g in genomes:

            # find all the alignments of u and v on genome g
            alns_on_g = g_to_aln_idx[g]
            u_alns = node_to_aln_idx[u] & alns_on_g
            v_alns = node_to_aln_idx[v] & alns_on_g

            # if we don't have alignments on both u, v on the genome, then continue
            if not (u_alns and v_alns):
                continue

            # there needs to be at least one position on g, for which, u,v are in close proximity
            # for the edge to be explained by genome g
            aln_found = True
            if check_proximity(u_alns, v_alns, alignments, assembly_artifact_gap):
                valid_edge = True
                break  # job done, edge can't be a fp

        # then unexplained by any genome
        if not valid_edge:
            print(asm_name, 'FPT1' if aln_found else 'FPT2', u, v)
            fp_record['value'] += 1


    # iterate through genomes and compute FN
    print('Computing TP...')
    for g in tqdm.tqdm(genomes):
        alns_on_g = [alignments[a] for a in g_to_aln_idx[g]]
        consensus_breakpoints = breakpoint_dict[g]
        ivl_tree = make_interval_tree(alns_on_g)
        for bp in consensus_breakpoints:
            left_bound, right_bound = bp - window_sz, bp + window_sz
            exists, path = path_exists(
                ivl_tree[left_bound: right_bound], left_bound, right_bound,
                asm.adjM, assembly_artifact_gap, record_path=False
            )
            if exists:
                tp_records[g]['value'] += 1
            else:
                fn_records[g]['value'] += 1

    total_records = list(tp_records.values()) + [fp_record] + list(fn_records.values())
    return total_records


def make_interval_tree(alns):
    intervals = []
    for node in alns:
        if node.start == node.end:
            continue
        intervals.append(Interval(node.start, node.end, node.contig_name))
    interval_tree = IntervalTree(intervals)
    return interval_tree


def path_exists(nodes, left_bound, right_bound, adjM, aag, record_path=True):
    assert record_path in [True, False]
    exists = False
    nodes = set(nodes)
    path = list() if record_path else None

    failed_from = set()
    for u in nodes:
        if u.begin > left_bound:
            continue
        exists = search(u, nodes - {u}, right_bound, adjM, failed_from, aag, path)
        if exists:
            return (exists, path) if record_path else (exists, None)
    return (exists, path) if record_path else (exists, None)


def undersampled_region(nodes, left_bound, right_bound):
    left_spanned, right_spanned = False, False
    for n in nodes:
        if n.begin <= left_bound:
            left_spanned = True
        if n.end >= right_bound:
            right_spanned = True

    return not (left_spanned and right_spanned)


def search(u, nodes, right_bound, adjM, failed_from, aag, path=None):
    # We already tried to search from this node and it failed,
    # so save time and avoid re-computing the same outcome
    if u in failed_from:
        return False

    if u.end >= right_bound:
        if path is not None:
            path.append((u, None))
        return True

    exists = False
    for v in nodes:
        if abs(v.begin - u.end) <= aag and v.end > u.end and (adjM.loc[u.data, v.data] or adjM.loc[v.data, u.data]):
            if path is not None:
                path.append((u, v))
            exists = search(v, nodes - {v}, right_bound, adjM, failed_from, aag, path)
        if exists:
            return exists

    # Returning from here means the right_bound was unreachable from
    # u, so, we failed to find the assembly from u. Record this in a list,
    # so other search call that extend through u will know its a waste of time
    # and terminate early
    failed_from.add(u)
    return exists


def check_proximity(u_alns, v_alns, alignments, aag):
    gap_distances = list()
    for u_aln, v_aln in product(u_alns, v_alns):
        u_aln, v_aln = alignments[u_aln], alignments[v_aln]
        assert(u_aln.start <= u_aln.end)
        assert(v_aln.start <= v_aln.end)
        # We don't know whether u comes before v or v before u in this pair.
        # Do we take the minimum of "gap" between end and start in both orientations
        gap_distance = min(abs(u_aln.end - v_aln.start), abs(u_aln.start - v_aln.end))
        gap_distances.append(gap_distance)
    #print(gap_distances)
    return min(gap_distances) <= aag


def compute_assembly_complexity(asms, key, dataset, n_pairs):
    # Count number of nodes and edges for each analysis and output to disk
    complexity = pd.DataFrame(index=asms.keys(), columns=['key', 'dataset', 'depth', 'num_nodes', 'num_edges'])
    for a in asms.values():
        complexity.loc[a.assembler, :] = [key, dataset, n_pairs, len(a.adjM.index), a.adjM.sum().sum()]
    return complexity


def good_alignment(hit, ctg_len, mode='multiple', min_similarity=0.98, length_epsilon=0.02):
    """
    Checks whether the alignment is good enough to be considered where the contig derived from in the genome.

    Within some epsilon, the length of the contig should be exactly the length of the alignment along the genome.
    The maximum allowed mismatches should also be (1 - min_similarity) % of the total sequence.
    """

    possible_modes = ['multiple']
    assert mode in possible_modes

    lower_bound, upper_bound = round((1 - length_epsilon) * ctg_len), round((1 + length_epsilon) * ctg_len)
    alignment_length = abs(hit.r_en - hit.r_st)
    if alignment_length < lower_bound or alignment_length > upper_bound:
        return False

    max_nm = int(ctg_len * (1 - min_similarity))
    if hit.NM <= max_nm:
        return True  # Then it's a good alignment

    return False


#def align_nodes(asm, aligner, g=None, all_alignments=True):
#    adjM, nodes = asm.adjM, asm.contigs
#    alignments = list()
#    # Construct the tour
#    unaligned_nodes = set()
#    for idx, (name, seq) in tqdm.tqdm(enumerate(nodes)):  # read a fasta/q sequenceA
#        no_hits = True
#        for hit in aligner.map(seq):  # traverse alignments
#            if g is not None and (hit.ctg != g):
#                continue
#            if all_alignments or good_alignment(hit, len(seq)):
#                no_hits = False
#                alignments.append(Alignment(name, hit.r_st, hit.r_en, hit.ctg, len(seq), hit.is_primary))
#        if no_hits:
#            unaligned_nodes.add(idx)
#    return alignments, unaligned_nodes

def align_nodes_para(aligner, g, all_alignments, dat):
    idx, (name, seq) = dat
    alignments = list()
    for hit in aligner.map(seq):  # traverse alignments
        if g is not None and (hit.ctg != g):
            continue
        if all_alignments or good_alignment(hit, len(seq)):
            alignments.append(Alignment(name, hit.r_st, hit.r_en, hit.ctg, len(seq), hit.is_primary))
    return alignments, idx

def align_nodes(asm, aligner, g=None, all_alignments=True):
    adjM, nodes = asm.adjM, asm.contigs
    alignments = list()
    unaligned_nodes = set()
    fun = functools.partial(align_nodes_para, aligner, g, all_alignments)
    with TPool(processes=16) as pool:
        ret = pool.map(fun, enumerate(nodes))
        for a, i in ret:
            if len(a) == 0:
                unaligned_nodes.add(i)
            else:
                alignments += a
    return alignments, unaligned_nodes
    


def compute_consensus_breakpoints(
        alignment_dict,
        genomes,
        genome_lens,
        max_within_clust_dist,
        window_sz,
        average=np.median,
        contigs_only=True,
        cluster=True
):
    """
    Constructs consensus breakpoints from multiple assemblies alignments.

    Each assembler has a set of alignments - intervals describing the alignment of it's
    nodes against a genome. The start/endpoints of these alignments are all flattened
    into a single list, breakpoints , which then clustered and an average point (consensus
    breakpoint) is computed per cluster.

    A consensus breakpoint is the average position of a set of breakpoints. It represents
    the consensus coordinate of a difficult to assemble region of a genome.
    """
    breakpoints_dict = defaultdict(list)
    for i, g in enumerate(genomes):
        covered_region = np.full(genome_lens[i], 0, dtype=np.uint32)
        breakpoints = list()
        for k, (alignments, u) in alignment_dict.items():
            if contigs_only and (k != 'megahit_contigs' and k != 'metaspades_contigs'):
                continue
            for a in alignments:
                if a.genome != g:
                    continue
                breakpoints.append(a.start)
                breakpoints.append(a.end)
                covered_region[a.start: a.end] += 1

        breakpoints = np.array(breakpoints)
        remaining_bps = np.full(len(breakpoints), True)
        for j, bp in enumerate(breakpoints):
            if (bp - window_sz < 0) or (bp + window_sz >= len(covered_region)):
                remaining_bps[j] = False
                continue
            if not ((covered_region[bp-window_sz] == 1) and (covered_region[bp+window_sz] == 1)):
                remaining_bps[j] = False
        breakpoints = breakpoints[remaining_bps]
        #print('breakpoints removed due to undersampling', len(remaining_bps) - sum(remaining_bps))
        #print('Total breakpoints for genome', g, sum(remaining_bps))
        if len(breakpoints) == 0:
            continue
        if cluster and len(breakpoints) > 1:
            breakpoints = np.array(list(breakpoints)).reshape(-1, 1)  # reshape needed to work with sk-learn expected input
            cluster_labels = AgglomerativeClustering(
                distance_threshold=max_within_clust_dist, n_clusters=None, linkage='complete'
            ).fit(breakpoints).labels_

            # compute consensus_breakpoints
            breakpoints = breakpoints.reshape(breakpoints.shape[0])  # return to series
            consensus_breakpoints = [
                int(average(breakpoints[cluster_labels == label])) for label in set(cluster_labels)
            ]
            breakpoints_dict[g] = sorted(consensus_breakpoints)
        else:
            breakpoints_dict[g] = sorted(list(breakpoints))
    return breakpoints_dict

# class EvaluateAssemblies:
#
#    def __init__(self, alignments, asms, breakpoints, assembly_artifact_gap, window_size):
#        self.alignments = alignments
#        self.asms = asms
#        self.breakpoints = breakpoints
#        self.window_size = window_size
#        self.assembly_artifact_gap = assembly_artifact_gap
#
#    def compute_per_genome_cov_metrics(self, genome_len, assembler, mode='collapse_within_node'):
#        """
#        Counts the coverage tp, fp, fn for a given assembly.
#
#        A basepair of the genome is assigned:
#            - tp if the basepair is covered by a "good alignment" exactly once. Note, should repeats occur
#              we can be in the situation where each repeat has multiple good alignments to it, or no repeat has
#              any good alignments to it. What justifies calling  the "exactly once" you ask, particularly in the case
#              of repeats?
#
#              Well, if the repeats are identical, the assembly
#              should collapse the sampled data into a single sequence, hence there would be exactly 1 contig for
#              all positions (and therefore no multiple alignments). If the repeats are not identical, the
#              assembly should keep them separate, mapping each to their own location with near 100% identity to their
#              own location (i.e we'd have one contig per repeat, each mapping to its own loci, and therefore no
#              multiple alignments). The gray area is due to phasing, where the assembly mixes non-identical basepairs
#              from both repeats. In this case, the contig(s) fail to map to any repeat location as a good alignment.
#            - fp if the basepair is covered more than once, so the assembly is spuriously repeating the genome.
#              This can happen with repeats, but, hopefully the above justified the calling a fp to you.
#            - fn if the basepair is covered 0 times, so the assembly is skipping parts of the genome, we have no data
#              captured here.
#        """
#        cov_vector = np.zeros(genome_len)
#        if mode == 'collapse_within_node':
#
#            # separate alignments by node name
#            alignments_by_node = defaultdict(list)
#            for a in self.alignments[assembler]:
#                alignments_by_node[a.name].append(a)
#
#            # Overlapping alignments within a node cannot increase the fp_score
#            # only overlapping alignments between nodes
#            for alns in alignments_by_node.values():
#                node_vector = np.zeros(genome_len)
#                for a in alns:
#                    node_vector[a.start: a.end] = 1
#                cov_vector += node_vector
#
#        else:
#            for node in self.get_unique_alignments(set(self.alignments[assembler])):
#                cov_vector[node.start: node.end] += 1
#
#        return np.sum(cov_vector == 1), np.sum(cov_vector > 1), np.sum(cov_vector == 0)
#
#    def get_unique_alignments(self, alignments, epsilon=100):
#
#        # sort alignments by length, largest to smallest
#        alignments = sorted(alignments, key=lambda x: abs(x.start - x.end), reverse=True)
#
#        # Cluster the start and ends of alignments
#        starts = np.array([a.start for a in alignments]).reshape(-1, 1)
#        ends = np.array([a.end for a in alignments]).reshape(-1, 1)
#        start_clust = AgglomerativeClustering(
#            distance_threshold=epsilon, n_clusters=None, linkage='complete'
#        ).fit(starts).labels_
#        end_clust = AgglomerativeClustering(
#            distance_threshold=epsilon, n_clusters=None, linkage='complete'
#        ).fit(ends).labels_
#
#        # Get the index of remaining alignments
#        remaining = set(range(len(alignments)))
#
#        # get index to each alignment
#        pos = np.arange(len(alignments))
#
#        # iterate through alignments, from largest to smallest
#        # and remove alignments that are similar
#        for i, a in enumerate(alignments):
#
#            # if a previous alignment already deleted the current, then
#            # skip over the current
#            if i not in remaining:
#                continue
#
#            # get the positions (indices) of alignments that have the same cluster
#            # assignment as the ith alignment for start and end. Then take the intersection,
#            # which identifies alignments with similar starts and ends to the ith alignment
#            similar = set(pos[end_clust == end_clust[i]]).intersection(set(pos[start_clust == start_clust[i]]))
#
#            for s in similar - {i}:
#                remaining.remove(s)
#
#        return [alignments[i] for i in remaining]
#
#
#    def compute_per_genome_cnx_metrics(self, assembler):
#
#        interval_tree = self.make_interval_tree(self.alignments[assembler])
#        tp_count, fn_count, fp_count = 0, 0, 0
#        adjM = self.asms[assembler].adjM
#        for bp in tqdm.tqdm(self.breakpoints):
#            left_bound, right_bound = bp - self.window_size, bp + self.window_size
#            exists, path = self.path_exists(
#                interval_tree[left_bound: right_bound], left_bound, right_bound, adjM, record_path=False
#            )
#
#            if exists:
#                tp_count += 1
#            else:
#                fn_count += 1
#
#        # Compute the false negatives
#        fp_count = self.compute_per_genome_cnx_fp(assembler, adjM)
#
#        return tp_count, fp_count, fn_count
#
#    def compute_per_genome_cnx_fp(self, assembler, adjM):
#
#        alignments = self.alignments[assembler]
#        r, c = np.where(np.tril(adjM.values))
#        aligned_nodes = {aln.name for aln in alignments}
#
#        # Get lookup table from node -> alignments
#        node_to_alignments = defaultdict(list)
#        for idx, a in enumerate(alignments):
#            node_to_alignments[a.name].append(idx)
#
#        # iterate through each edge and determine whether the should exist based on node
#        # proximity
#        fp_count_aligned, fp_count_inf = 0, 0
#        for u, v in zip(adjM.index[r], adjM.index[c]):
#
#            u_aligned, v_aligned = u in aligned_nodes, v in aligned_nodes
#            # Both contigs mapped to the genome
#            if u_aligned and v_aligned:
#
#                # Since u and v could align to multiple positions in the genome.
#                # If the distance between the closest instances of u and v is too
#                # large, then we have a false positive
#                u_alns = node_to_alignments[u]
#                v_alns = node_to_alignments[v]
#                gap_distances = list()
#                for u_aln, v_aln in product(u_alns, v_alns):
#                    u_aln, v_aln = alignments[u_aln], alignments[v_aln]
#                    # We don't know whether u comes before v or v before u in this pair.
#                    # Do we take the minimum of "gap" between end and start in both orientations
#                    gap_distance = min(abs(u_aln.end - v_aln.start), abs(u_aln.start - v_aln.end))
#                    gap_distances.append(gap_distance)
#                if min(gap_distances) > self.assembly_artifact_gap:
#                    # print('FPEDGE', u, v)
#                    fp_count_aligned += 1
#
#            # Only one contig aligned to genome, assume minimum gap distance between u and v is inf
#        #            elif (u_aligned and not v_aligned) or (not u_aligned and v_aligned):
#        #                print('FPEDGEINF', u, v)
#        #                fp_count_inf += 1
#
#        print(assembler, 'fp_count_aligned', fp_count_aligned)
#        return fp_count_aligned
#
# def compute_cnx_fp_from_unaligned_contigs(aligner, contigs, adjM, min_len=300):
#    """
#    For contigs u, v that fail to align to any genome. An edge between such contigs is considered a false positive.
#    """
#
#    # Collect all contigs that failed to align
#    fp = 0
#    unaligned = set()
#    for name, seq in contigs:
#        if len(seq) < min_len:
#            continue
#        hits = aligner.map(seq)
#        for hit in hits:
#            if not good_alignment(hit, len(seq)):
#                unaligned.add(name)
#
#    for i, j in combinations(unaligned, 2):
#        if adjM.loc[i, j] or adjM.loc[j, i]:
#            fp += 1
#    return fp


# def compute_cnx_fp_from_inter_genome_edges(aligner, asm, alignments_by_genome, min_len=300):
#    """
#    An edge between two contigs u, v is only permitted if u, v align to at least one common genome.
#    Otherwise, the does not exist a genome supporting the edge.
#
#    """
#
#    # get location of all edges
#    r, c = np.where(np.tril(asm.adjM.values))
#
#    # get all alignments as a lookup from nodes -> alignment
#    aligned_nodes = {g: alignments_by_genome[g][asm.assembler] for g in aligner.seq_names}
#    fp = 0
#    for u, v in zip(asm.adjM.index[r], asm.adjM.index[c]):
#        explained = False
#        for g in aligner.seq_names:
#            aligned_to_g = {a.name for a in aligned_nodes[g]}
#            u_aligned, v_aligned = u in aligned_to_g, v in aligned_to_g
#            # Both contigs mapped to the genome. If they map to the genome, but not close enough
#            # they will have been counted as a per-genome fp already, so we are not checking that case here.
#            if u_aligned and v_aligned:
#                explained = True
#
#        # then there does not exist a genome that u, v both align to. This means the edge can only be explaind
#        # by a chimeric event, and is therefore a false positive
#        if not explained:
#            fp += 1
#    print('fp inter genome edge', fp)
#    return fp
# def compute_cov_fp_from_unaligned_contigs(aligner, contigs):
#    """
#        We consider all contigs that fail to have a good alignment to any position in any of the genomes as
#        false positive contig sequences. These are simply too noisy to align to any position, so they fail
#        to be a faithful representation of any position in the genome.
#
#    """
#
#    # We also consider all the contigs that have not good alignment as fp_counts
#    fp_count = 0
#    for name, seq in contigs:
#        # map the contig against all references
#        hits = aligner.map(seq)
#        ga = False
#        # check to see if at least one of the hits is a good alignment, and if it is, break
#        for hit in hits:
#            if good_alignment(hit, len(seq)):
#                ga = True
#                break
#        if not ga:
#            fp_count += len(seq)
#    return fp_count




#def get_coverage_metrics(asm_name, alignments, genomes, genome_lens, asm, node_to_aln_index, g_to_aln_index, unassigned_contigs, minimap_kmer_size=15):
#    # get record objects
#    tp_records = {
#        g: {'genome': g, 'assembler': asm_name, 'metric': 'cov_tp', 'value': 0} for g in genomes
#    }
#
#    # we can have a per-genome fp, as sequences that "overcover" a genome, as well as non-genome-specific fp, which
#    # fail to map anywhere
#    fp_records = {
#        g: {'genome': g, 'assembler': asm_name, 'metric': 'cov_fp', 'value': 0} for g in genomes
#    }
#    fp_records['-'] = {'genome': '-', 'assembler': asm_name, 'metric': 'cov_fp', 'value': 0}
#
#    fn_records = {
#        g: {'genome': g, 'assembler': asm_name, 'metric': 'cov_fn', 'value': 0} for g in genomes
#    }
#
#    for i, g in tqdm.tqdm(enumerate(genomes)):
#        # filter alignments for those of high quality. We want to measure the coverage of contigs against a genome
#        # by aligning contigs to the position they derive from, hence we pick the location of minimal edit distance for
#        # each contig
#        aln_on_g = g_to_aln_index[g]
#        if len(aln_on_g) == 0:
#            continue
#        cov_vector = np.zeros(genome_lens[i])
#        node_vector = np.zeros(genome_lens[i])
#        for node, alns in node_to_aln_index.items():
#            if len(alns) == 0:
#                continue
#            node_vector[:] = 0  # reset to 0
#            alns = alns & aln_on_g  # intersect with best alignments on the genome
#            for aln_idx in alns:
#                a = alignments[aln_idx]
#                node_vector[a.start: a.end] = 1
#            cov_vector += node_vector
#        tp_records[g]['value'] = np.sum(cov_vector > 0)
#        fn_records[g]['value'] = np.sum(cov_vector == 0)
#        fp_records[g]['value'] = np.sum(cov_vector > 1)
#
#    # contigs that fail to align to anny genome
#    num_unassigned_bp = 0
#    for idx in unassigned_contigs:
#        n, seq = asm.contigs[idx]
#        if len(seq) < minimap_kmer_size:
#            continue
#        num_unassigned_bp += len(seq)
#    fp_records['-']['value'] = num_unassigned_bp
#
#    return list(tp_records.values()) + list(fn_records.values()) + list(fp_records.values())
