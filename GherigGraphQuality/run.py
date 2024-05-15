
import yaml
import os
import sys
from collections import defaultdict
import pandas as pd
import mappy as mp

from Utils.AssemblyParser import Assembly
from Utils.evaluate_assembly import \
    align_nodes, compute_assembly_complexity, get_connectivity_metrics, \
    get_coverage_metrics, get_top_alignments, compute_consensus_breakpoints, compute_assembly_nX

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: <exe> <config>')
        sys.exit()

    # parse configs
    with open(sys.argv[1]) as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        asm_artifact_gap = int(config['asm_artifact_gap'])
        max_within_clust_distance = int(config['max_within_clust_distance'])
        window_size = config['window_size']
        reference = config['reference']
        out_dir = config['out_dir']
        run_desc = config['run_description']
        key = config['key']
        dataset = config['dataset']
        depth = -1
        
    # make out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for k, v in config.items():
        print(k, v)

    # Parse each tool's assembly file output into a common format of an adjacency matrix and assembly sequences
    for k, v in config['ASMS'].items():
        print(k, v)
    asms = [Assembly(**v) for v in config['ASMS'].values()]
    asms = {a.assembler: a for a in asms}

    # Compute complexity and nX
    #complexity = compute_assembly_complexity(asms, key, dataset, depth)
    #complexity.to_csv(os.path.join(out_dir, f'{run_desc}_complexity.csv'))
    #nX = compute_assembly_nX(asms, key, dataset, depth)
    #nX.to_csv(os.path.join(out_dir, f'{run_desc}_nX.csv'))
    #print(complexity)
    #print(nX)

    # Build a minimap alignment object indexing the metagenomic reference genomes
    aligner = mp.Aligner(reference)

    # get genome information from alignment tool
    genomes = list(aligner.seq_names)
    genome_lens = [len(aligner.seq(g)) for g in genomes]

    # compute alignments and breakpoints
    alignments_dict = dict()
    for asm_name, asm in asms.items():
        print('Computing alignments for assembly', asm_name, '...')
        alignments, unaligned_contigs = align_nodes(asm, aligner, all_alignments=False)
        alignments_dict[asm_name] = (alignments, unaligned_contigs)

    print('Computing breakpoints...')
    breakpoint_dict = compute_consensus_breakpoints(
        alignments_dict,
        genomes,
        genome_lens,
        max_within_clust_distance,
        window_size
    )

    # Run evaluation per assembler
    stats_dfs = list()
    for asm_name, asm in asms.items():
        print(f'Evaluating assembler {asm_name}')

        # get all alignments for assembler
        alignments, unaligned_contig_indicies = alignments_dict[asm_name]

        # build genome to alignment lookup
        g_to_aln_idx = defaultdict(set)
        for i, aln in enumerate(alignments):
            g_to_aln_idx[aln.genome].add(i)

        # build node to alignment lookup
        node_to_aln_idx = defaultdict(set)
        for i, aln in enumerate(alignments):
            node_to_aln_idx[aln.contig_name].add(i)


        # get connectivity metrics
        print('Computing connectivity metrics')
        cnx_metric_records = get_connectivity_metrics(
            asm_name, alignments, genomes, g_to_aln_idx, node_to_aln_idx, breakpoint_dict,
            asm, asm_artifact_gap,  window_size,
            unaligned_contig_indicies
        )

        # Filter the alignments for the best hit per contig. Intersect all per-genome alignments with the top alignments
        # Equally good top alignments are retained
        top_alignment_per_contig = get_top_alignments(alignments)
        for g in g_to_aln_idx:
            g_to_aln_idx[g] = g_to_aln_idx[g] & top_alignment_per_contig

        # build contig map
        contig_map = None
        print('Computing coverage metrics')
        cov_metric_records = get_coverage_metrics(
            asm_name, alignments, genomes, genome_lens, asm, g_to_aln_idx, unaligned_contig_indicies, contig_map
        )

        # combine records and add additional metadata
        all_records = cnx_metric_records + cov_metric_records
        all_records = [{**e, 'key':key, 'dataset':dataset, 'depth': depth} for e in all_records]

        # convert to table
        df = pd.DataFrame().from_records(all_records)
        stats_dfs.append(df)

    # write out data
    stats_dfs = pd.concat(stats_dfs)
    stats_dfs.to_csv(os.path.join(out_dir, f'{run_desc}_quality.csv'), index=False)
