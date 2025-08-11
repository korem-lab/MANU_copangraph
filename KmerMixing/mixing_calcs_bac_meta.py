
import os
import sys
import mappy as mp
import pandas as pd
from glob import glob
from Utils.evaluate_assembly import get_coverage_metrics, align_nodes, n50, num_edges, num_nodes, get_top_alignments
from Utils.AssemblyParser import Assembly
from collections import defaultdict
import re

#COPANGRAPHS = './data/KmerMixing/gherig/copangraph/*.gfa'
#MEGAHIT_GRAPHS = './data/KmerMixing/gherig/coassembly/megahit/intermediate_contigs/*.fastg'

OUTDIR = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta'

def parse_name(s, gherig=True):
    if gherig:
        if 'fastg' in s:
            return 'megahit', float(re.findall('k([0-9]+)',s)[0])
            
        else:
            return 'copangraph', float(re.findall('(0.[0-9]+)',s)[0])
    else:
        if 'fastg' in s:
            fields = s.split('/')
            replicate, gen, rest = fields[-3].split('_')
            return 'megahit', replicate, int(gen), int(fields[-1][1:].replace('.fastg', ''))
        else:
            replicate, gen, sd, rest = os.path.basename(s).split('_')
            return 'copangraph', replicate, int(gen), float(sd)

def total_aligned_bp(df):
    df = df.sort_values(['genome','start']).reset_index(drop=True)
    group = df.groupby('genome').apply(lambda x: (x['start'] > x['end'].cummax().shift(fill_value=-float('inf'))).cumsum(), include_groups=False).reset_index()
    print(group)
    ivl_merge = df.groupby([group.genome, group[0]]).agg({'start':'min', 'end':'max'}).reset_index()
    return (ivl_merge.end - ivl_merge.start).sum()

if __name__ == '__main__':

    graph, window_sz, genomes = sys.argv[1], sys.argv[2], sys.argv[3]
    window_sz = int(window_sz)
    aligner = mp.Aligner(genomes, k=15, w=window_sz)
    genomes = list(aligner.seq_names)
    genome_lens = [len(aligner.seq(g)) for g in genomes]
    cov_records = list()
    mgs_records = list()

    name = os.path.basename(graph.replace('.fasta', '').replace('.gfa', ''))
    asm = Assembly(assembler=name, assembly_file=graph)
    alignments, unaligned_contigs = align_nodes(asm, aligner, mag=True)
    g_to_aln_idx = defaultdict(set)
    for i, aln in enumerate(alignments):
        g_to_aln_idx[aln.genome].add(i)
    contig_map = {name: i for i, name in enumerate(set(n for n,_ in asm.contigs))}
    cov_record = get_coverage_metrics(
        name, alignments, genomes, genome_lens, asm, g_to_aln_idx, unaligned_contigs, contig_map, disregard_unmapped=False
    )
    cov_record = pd.DataFrame(cov_record)
    cov_record = cov_record.groupby('metric').value.sum()
    prc = cov_record.cov_tp / (cov_record.cov_tp + cov_record.cov_fp)
    rcl = cov_record.cov_tp / (cov_record.cov_tp + cov_record.cov_fn)
    fsc = (2*prc*rcl) / (prc+rcl)
    n_50 = n50(graph)
    nodes = num_nodes(asm)
    edges = num_edges(asm)
    tool, rep, gen, param= parse_name(graph, gherig=False) 
    cov_records.append([rep, gen, tool, param, name, n_50, nodes, edges, 'cov_tp', cov_record.cov_tp])
    cov_records.append([rep, gen, tool, param, name, n_50, nodes, edges, 'cov_fn', cov_record.cov_fn])
    cov_records.append([rep, gen, tool, param, name, n_50, nodes, edges, 'cov_fp', cov_record.cov_fp])
    cov_records.append([rep, gen, tool, param, name, n_50, nodes, edges, 'cov_precsion', prc])
    cov_records.append([rep, gen, tool, param, name, n_50, nodes, edges, 'cov_recall', rcl])
    cov_records.append([rep, gen, tool, param, name, n_50, nodes, edges, 'cov_F-score', fsc])
    cov_records = pd.DataFrame(cov_records, columns=['rep', 'gen', 'assembler', 'parameter', 'name', 'n50', 'num_nodes', 'num_edges', 'metric', 'value'])

    # convert named tuples to df
    aln_df = pd.DataFrame(alignments)

    is_multi_genome = aln_df.groupby('contig_name').apply(lambda x: x.genome.nunique()>1)
    aln_df_multigenome  = aln_df.loc[aln_df.contig_name.isin(is_multi_genome.index[is_multi_genome]),:]
    if is_multi_genome.sum() != 0:
        multigenome_coverage_bp = total_aligned_bp(aln_df_multigenome)
    else:
        multigenome_coverage_bp = 0
    max_len = asm.get_max_node_lens()
    max_len = max_len.reset_index()
    max_len.columns = ['node_id', 'max_len']
    sum_maxlen_bp = max_len.max_len.sum()

    # reduce to contigs
    multi_genome_status = aln_df.groupby('contig_name').apply(
        lambda x: pd.Series(
            [rep, gen, tool, param, name, n_50, nodes, edges, 
             x.genome.nunique() > 1, x.shape[0], multigenome_coverage_bp, sum(genome_lens), sum_maxlen_bp], 
            index=['rep', 'gen', 'assembler', 'parameter', 'name', 'n50', 'num_nodes', 'num_edges', 
                   'is_multi_genome', 'num_alignments', 'multigenome_coverage_bp', 'total_ref_bps', 'sum_maxlen_total_bps']
        ),
        include_groups=False
    )
    multi_genome_status = multi_genome_status.reset_index(names='node_id')    
    multi_genome_status = pd.merge(multi_genome_status, max_len, on='node_id', how='inner')


    mgs_record = pd.DataFrame([[
        rep, gen, tool, param, name, n_50, nodes, edges, sum_maxlen_bp, sum(genome_lens),
        multi_genome_status.is_multi_genome.sum(),  # proportion nodes multigenome
        multi_genome_status.is_multi_genome.sum()/nodes,  # proportion nodes multigenome
        multi_genome_status.loc[multi_genome_status.is_multi_genome, 'max_len'].sum(), # sum_maxlen_multigenome
        multi_genome_status.loc[multi_genome_status.is_multi_genome, 'max_len'].sum()/sum_maxlen_bp, # proportion maxlen multigenome
        multigenome_coverage_bp, # multigenome coverage
        multigenome_coverage_bp/sum(genome_lens) # proportion multigenome coverage of reference
    ]], columns=[
        'rep', 'gen', 'assembler', 'param', 'name', 'n50', 'num_nodes', 'num_edges', 'sum_maxlen_bp', 'total_ref_bps',
        'sum_multigenome_nodes', 'prop_multigenome_nodes', 'sum_maxlen_multigenome_bp', 'prop_maxlen_multigenome_bp', 'sum_multigenome_coverage_bp', 'prop_multigenome_coverage_bp'
    ])

    cov_records.to_csv(os.path.join(OUTDIR, f'{rep}_{gen}_{tool}_{param}_coverage_records.csv'), index=False)
    print(mgs_record)
    mgs_record.to_csv(os.path.join(OUTDIR, f'{rep}_{gen}_{tool}_{param}_multisample_record.csv'), index=False)