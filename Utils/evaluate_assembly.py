import numpy as np
import Utils.parse_seq as ps
import pandas as pd
from intervaltree import IntervalTree, Interval
from collections import defaultdict, namedtuple
from itertools import product, combinations
from sklearn.cluster import AgglomerativeClustering
from scipy import sparse as sp
import tqdm
from multiprocessing.dummy import Pool as TPool
from multiprocessing import Pool as PPool
from concurrent.futures import ThreadPoolExecutor as TP
import functools

Alignment = namedtuple('Alignment', ['contig_name', 'start', 'end', 'genome', 'contig_len', 'is_primary', 'sam'])

def rev_comp(s):
    return ''.join(
        {'A':'T', 'T':'A', 'C':'G','G':'C'}[e] for e in s[::-1]
    )

def paf_to_sam(paf_line, seq):
    """
    Convert a single PAF line to a simplified SAM line.

    :param paf_line: str, a single line from a PAF file
    :return: str, corresponding SAM line
    """
    fields = paf_line.strip().split('\t')
    if len(fields) < 12:
        raise ValueError("PAF line must have at least 12 fields")

    qname = fields[0].replace(' ', '_').replace('flag=', 'f')
    qlen = int(fields[1])
    qstart = int(fields[2])+1
    qend = int(fields[3])+1
    strand = fields[4]
    rname = fields[5]
    rlen = int(fields[6])
    rstart = int(fields[7])+1
    rend = int(fields[8])+1
    match_len = int(fields[9])
    aln_len = int(fields[10])
    mapq = int(fields[11])
    cigar = fields[14].split(':')[-1]


    # Strand
    flag = 0 if strand == '+' else 16

    # POS: SAM is 1-based, PAF is 0-based
    pos = rstart

    # CIGAR: simplified to match_len + softclips (not exact)
    # This is highly simplified and not always accurate

    # RNEXT, PNEXT, TLEN are unknown
    rnext = "*"
    pnext = 0
    tlen = qlen

    # Sequence and quality are unknown in PAF
    qual = "*"

    sam_fields = [
        qname,         # QNAME
        str(flag),     # FLAG
        rname,         # RNAME
        str(pos),      # POS
        str(mapq),     # MAPQ
        cigar,         # CIGAR
        rnext,         # RNEXT
        str(pnext),    # PNEXT
        str(tlen),     # TLEN
        seq[qstart-1:qend-1] if strand == '+' else rev_comp(seq[qstart-1:qend-1]),           # SEQ
        qual           # QUAL
    ]

    return '\t'.join(sam_fields)

def get_top_alignments(alignments, mode='is_primary'):
    assert mode in [
        'is_primary', 'mapq', 'AS', 'accept_all'
    ]

    if mode == 'accept_all':
        return set(range(len(alignments)))
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


def get_coverage_metrics(asm_name, alignments, genomes, genome_lens, asm, g_to_aln_index, unassigned_contigs, contig_map, minimap_kmer_size=15, disregard_unmapped=True):
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

    #fun = functools.partial(get_genome_cov, g_to_aln_index, genome_lens, contig_map, alignments) 
    total_bp = 0
    covered_bp = 0
    for (i, g) in tqdm.tqdm(enumerate(genomes)):
        g, tp, fn, fp = get_genome_cov(g_to_aln_index, genome_lens, contig_map, alignments, (i, g))
        covered_bp += tp
        total_bp += genome_lens[i]
        tp_records[g]['value'] = tp
        fn_records[g]['value'] = fn
        fp_records[g]['value'] = fp
    print('PROP COVERED: ', covered_bp/total_bp)

    # contigs that fail to align to any genome
    if not disregard_unmapped:
        num_unassigned_bp = 0
        for idx in unassigned_contigs:
            n, seq = asm.contigs[idx]
            if len(seq) < minimap_kmer_size:
                continue
            num_unassigned_bp += len(seq)
        fp_records['-']['value'] = num_unassigned_bp

    return list(tp_records.values()) + list(fn_records.values()) + list(fp_records.values())
        
def get_genome_cov(g_to_aln_index, genome_lens, contig_map, alignments, dat):
    g_idx, g = dat
    aln_on_g = sorted(list(g_to_aln_index[g]))
    if len(aln_on_g) == 0:
        return (g, 0, genome_lens[g_idx], 0)
    if contig_map is None:
        cov_vector = np.zeros(genome_lens[g_idx], dtype=np.uint32)
        for i, aln_idx in enumerate(aln_on_g):
            a = alignments[aln_idx]
            cov_vector[a.start:a.end] += 1
    else:
        cov_vector = np.zeros(genome_lens[g_idx], dtype=np.uint32)
        aln_on_g = sorted([alignments[i] for i in aln_on_g], key=lambda x: x.contig_name)
        i = 0
        if len(aln_on_g) == 1:
            a = alignments[i]
            cov_vector[a.start:a.end] += 1
        else:
            while i < len(aln_on_g)-1:
                temp_vector = np.zeros(genome_lens[g_idx], dtype=np.uint32)
                temp_vector[aln_on_g[i].start: aln_on_g[i].end] = 1
                while i + 1 < len(aln_on_g) and aln_on_g[i].contig_name == aln_on_g[i+1].contig_name:
                    temp_vector[aln_on_g[i+1].start: aln_on_g[i+1].end] = 1
                    i += 1
                cov_vector = cov_vector + temp_vector
                temp_vector = np.zeros(genome_lens[g_idx], dtype=np.uint32)
                i += 1
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
        unaligned,
        mh_graph_asm=None,
        mh_contig_to_node=None,
        copan_node_to_contig=None,
        mode='never_disregard_unmapped'
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
    total_edges = 0
    both_unmapped = 0
    one_unmapped = 0
    both_mapped =0


    def edge_in_megahit_graph(name, u, v, mh_graph_asm, copan_node_to_contig, mh_contig_to_node):
        if mh_graph_asm is None:
            return True
        if name == 'copangraph':
            u_contig = copan_node_to_contig[u]
            v_contig = copan_node_to_contig[v]
            u_mh = mh_contig_to_node[u_contig]
            v_mh = mh_contig_to_node[v_contig]
            u_mh = mh_graph_asm.adjM.index[mh_graph_asm.adjM.index.to_series().apply(lambda x: f'{u_mh[1:]}_' in x)][0]
            v_mh = mh_graph_asm.adjM.index[mh_graph_asm.adjM.index.to_series().apply(lambda x: f'{v_mh[1:]}_' in x)][0]
            return mh_graph_asm.adjM.loc[u_mh, v_mh] or mh_graph_asm.adjM.loc[v_mh, u_mh]
        elif name == 'megahit_graph':
            return mh_graph_asm.adjM.loc[u, v] or mh_graph_asm.adjM.loc[v,u]
        else:
            return False


    for u, v in tqdm.tqdm(zip(asm.adjM.index[r], asm.adjM.index[c])):
        total_edges += 1

        # If the edge involves an unaligned contig can
        # chose to consider it a false positive by default, or consider it unmeasureable
        u_aln = node_to_aln_idx[u]
        v_aln = node_to_aln_idx[v]


        if len(v_aln) == 0 and len(u_aln) == 0:
            both_unmapped += 1
        elif (len(v_aln) != len(u_aln)) and (len(u_aln) == 0 or len(v_aln) == 0):
            one_unmapped += 1
        else:
            both_mapped += 1

        if mode == 'disregard_if_any_unmapped':
            # If node u or node v has no mapped sequences
            # then cannot measure if fp status
            if len(v_aln) == 0 or len(u_aln) == 0:
                continue
        elif mode == 'disregard_if_both_unmapped':
            # if both u and v are unmapped, disregard, otherwise,
            # call it a false positve if one is mapped
            if len(v_aln) == 0 and len(u_aln) == 0:
                continue
            elif len(v_aln) == 0 ^ len(u_aln) == 0:
                fp_record['value'] += 1
        elif mode == 'never_disregard_unmapped':
            if len(v_aln) == 0 or len(u_aln) == 0:
                fp_record['value'] += 1
                if not edge_in_megahit_graph(asm_name, u, v, mh_graph_asm, copan_node_to_contig, mh_contig_to_node):
                    print('FP NOT IN MH UNMAPPED: ', u, v)
                continue
        

        valid_edge = False
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
            contiguous, pos = check_proximity(u_alns, v_alns, alignments, assembly_artifact_gap)
            if contiguous:
                valid_edge = True
                break  # job done, edge can't be a fp

        # then unexplained by any genome
        if not valid_edge:
            fp_record['value'] += 1
            if not edge_in_megahit_graph(asm_name, u, v, mh_graph_asm, copan_node_to_contig, mh_contig_to_node):
                print('FP NOT IN MH MAPPED: ', u, v)
    if total_edges != 0:
        print('Proportion both mapped: ', both_mapped/total_edges)
        print('Proportion neither mapped: ', both_unmapped/total_edges)
        print('Proportion one mapped: ', one_unmapped/total_edges)
    # iterate through genomes and compute TP/FN
    print('Computing TP...')
    tp_bed_records = list()
    fn_bed_records = list()
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
                tp_bed_records.append(f'{g}\t{bp}\t{bp+1}')
            else:
                fn_records[g]['value'] += 1
                fn_bed_records.append(f'{g}\t{bp}\t{bp+1}')

    total_records = list(tp_records.values()) + [fp_record] + list(fn_records.values())
    return total_records, tp_bed_records, fn_bed_records


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


def check_proximity(u_alns, v_alns, alignments, aag, asm_name=None):
    gap_distances = list()
    for u_aln, v_aln in product(u_alns, v_alns):
        u_aln, v_aln = alignments[u_aln], alignments[v_aln]
        assert(u_aln.start <= u_aln.end)
        assert(v_aln.start <= v_aln.end)
        # We don't know whether u comes before v or v before u in this pair.
        # Do we take the minimum of "gap" between end and start in both orientations
        #gap_distance = min(abs(u_aln.end - v_aln.start), abs(u_aln.start - v_aln.end))
        gd = abs(u_aln.end - v_aln.start)
        pos = (u_aln.end, v_aln.start)
        if abs(u_aln.start - v_aln.end) < gd:
            gd = abs(u_aln.start - v_aln.end) 
            pos = (u_aln.start, v_aln.end)
        gap_distances.append((gd, pos))
    if asm_name:
        print(f'{asm_name} GAPDIST:', min(gap_distances))
    #return min(gap_distances) <= aag
    mgd, mpos = sorted(gap_distances, key=lambda x: x[0])[0]
    return mgd <= aag, mpos


def num_nodes(a):
    return len(a.adjM.index)

def num_edges(a):
    return sp.tril(a.adjM.sparse.to_coo()).nnz


def good_alignment(hit, ctg_len, min_similarity=0.98, length_epsilon=0.02, mag=False):
    """
    Checks whether the alignment is good enough to be considered where the contig derived from in the genome.

    Within some epsilon, the length of the contig should be exactly the length of the alignment along the genome.
    The maximum allowed mismatches should also be (1 - min_similarity) % of the total sequence.
    """


    lower_bound, upper_bound = round((1 - length_epsilon) * ctg_len), round((1 + length_epsilon) * ctg_len)
    alignment_length = abs(hit.r_en - hit.r_st)
    if not (lower_bound <= alignment_length <= upper_bound):
        return False
    if mag:
        min_similarity = 0.999
    max_nm = int(ctg_len * (1.0 - min_similarity))
    if hit.NM <= max_nm:
        return True  # Then it's a good alignment

    return False

def align_nodes_para(aligner, g, mag, dat):
    idx, (name, seq) = dat
    alignments = list()
    for hit in aligner.map(seq):  # traverse alignmentsfor y in /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/yamls/*_analysis.yaml; do  sbatch --export=ALL,ANALYSIS=$y CAMISIMGraphQuality/batch_eval.sh ; done
        if g is not None and (hit.ctg != g):
            continue
        if good_alignment(hit, len(seq), mag=mag):
            alignments.append(Alignment(name, hit.r_st, hit.r_en, hit.ctg, len(seq), hit.is_primary, paf_to_sam(f'{name}\t{len(seq)}\t' + str(hit), seq)))
    return alignments, idx

def align_nodes(asm, aligner, g=None, mag=False):
    adjM, nodes = asm.adjM, asm.contigs
    alignments = list()
    unaligned_nodes = set()
    fun = functools.partial(align_nodes_para, aligner, g, mag)
    with TP(max_workers=4) as pool:
        ret = pool.map(fun, enumerate(nodes))
    for a, i in ret:
        if len(a) == 0:
            unaligned_nodes.add(i)
        else:
            alignments += a
    return alignments, unaligned_nodes

def write_bed(breakpoints, outfile):

    with open(outfile, 'w') as f:
        for g, bps in breakpoints.items():
            for bp in bps:
                f.write(f'{g}\t{bp}\t{bp+1}\n')

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
            if contigs_only and (k != 'megahit_contigs' and k != 'metaspades_contigs' and k != 'bcalm_contigs'):
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

def n50(asm):
    if type(asm) == str:
        if asm.endswith('.gfa'):
            is_copan = True
        else:
            is_copan = False
    
        with open(asm) as f:
            if is_copan:
                data = list(e for e in ps.parse_gfa(f) if e.type == ps.GFATypes.S)
            else:
                data = list(ps.parse(f, ps.Fasta))
        data = [e.seq for e in data]
    else:
        data = [e for _, e in asm.contigs]
    # organize sequences by node
    data = sorted(data, key=lambda x: len(x), reverse=True)
    total_50 = sum(len(e) for e in data)/2
    cum_sum= 0 
    for e in data:
        cum_sum += len(e)
        if cum_sum >= total_50:
            return len(e)
