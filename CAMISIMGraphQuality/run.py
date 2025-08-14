
import yaml
import os
import sys
from collections import defaultdict
import pandas as pd
import mappy as mp
from Utils.AssemblyParser import Assembly
from Utils.constants import depth_map
from Utils.evaluate_assembly import \
    align_nodes, get_connectivity_metrics, \
    get_coverage_metrics, get_top_alignments, compute_consensus_breakpoints

BURG = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/'
CPN = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/CPN'
MANITOU='/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/'
if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: <exe> <config>')
        sys.exit()

    # parse configs
    with open(sys.argv[1]) as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        config['reference'] = config['reference']
        config['out_dir'] = config['out_dir']
        asm_artifact_gap = int(config['asm_artifact_gap'])
        max_within_clust_distance = int(config['max_within_clust_distance'])
        window_size = config['window_size']
        reference = config['reference']
        out_dir = config['out_dir']
        key = config['key']
        run_desc = config['run_description']
        depth = depth_map(config['depth'][:-1])
        dataset, base_asm, _ = key.split('_')
        dataset = int(dataset)
        
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

    # Build a minimap alignment object indexing the metagenomic reference genomes
    aligner = mp.Aligner(reference)

    # get genome information from alignment tool
    genomes = list(aligner.seq_names)
    genome_lens = [len(aligner.seq(g)) for g in genomes]

    # compute alignments and breakpoints
    alignments_dict = dict()
    for asm_name, asm in asms.items():
        print('Computing alignments for assembly', asm_name, '...')
        alignments, unaligned_contigs = align_nodes(asm, aligner)
        unaligned_bp = sum(len(asm.contigs[i][1]) for i in unaligned_contigs)
        total_bp = sum(len(e[1]) for e in asm.contigs)
        print('TOTAL PROPORTION ALIGNED:', 1- unaligned_bp/total_bp)
        alignments_dict[asm_name] = (alignments, unaligned_contigs)

    for name, (alns,_) in alignments_dict.items():
        with open(os.path.join(out_dir, f'{run_desc}_{name}.alignment.sam'), 'w') as f:
            for g, gl in zip(genomes, genome_lens):
                f.write(f'@SQ\tSN:{g}\tLN:{gl}\n')
            for a in alns:
                f.write(a.sam + '\n')
            

    print('Computing breakpoints...')
    breakpoint_dict = compute_consensus_breakpoints(
        alignments_dict,
        genomes,
        genome_lens,
        max_within_clust_distance,
        window_size
    )
    #write_bed(breakpoint_dict, os.path.join(out_dir, f'{run_desc}_breakpoints.bed'))
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
        cnx_metric_records, tp_bed_records, fn_bed_records= get_connectivity_metrics(
            asm_name, alignments, genomes, g_to_aln_idx, node_to_aln_idx, breakpoint_dict,
            #asm, asm_artifact_gap,  window_size, # CURRENT MANUSCRIPT Jul 6
            asm, window_size*2,  window_size, # TESTING
            unaligned_contig_indicies, 
            mode='never_disregard_unmapped'
        )
        ## write tp and fn bed
        #with open(os.path.join(out_dir, f'{run_desc}_{asm_name}_tp.bed'), 'w') as f:
        #    f.write('\n'.join(tp_bed_records))
        #with open(os.path.join(out_dir, f'{run_desc}_{asm_name}_fn.bed'), 'w') as f:
        #    f.write('\n'.join(fn_bed_records))

        # Filter the alignments for the best hit per contig. Intersect all per-genome alignments with the top alignments
        # Equally good top alignments are retained
        top_alignment_per_contig = get_top_alignments(alignments, mode='accept_all')
        for g in g_to_aln_idx:
            g_to_aln_idx[g] = g_to_aln_idx[g] & top_alignment_per_contig

        # build contig map
        contig_map = {name: i for i, name in enumerate(set(n for n,_ in asm.contigs))}
        print('Computing coverage metrics')
        cov_metric_records = get_coverage_metrics(
            asm_name, alignments, genomes, genome_lens, asm, g_to_aln_idx, 
            unaligned_contig_indicies, contig_map, disregard_unmapped=False
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
