import os
import sys
import mappy as mp
import pandas as pd
from glob import glob
from ..Utils.evaluate_assembly import get_coverage_metrics, align_nodes, n50, num_edges, num_nodes
from ..Utils.AssemblyParser import Assembly
from collections import defaultdict
import re

GRAPHS='./graphdir'
GENOMES='./genomes'
OUTDIR = './OUT'

def parse_name(s):
    if 'fastg' in s:
        return 'megahit', float(re.findall('k([0-9]+)')[0])
    else:
        return 'copangraph', float(re.findall('sd_(0\.[0-9]+)')[0])

if __name__ == '__main__':

    graphs = glob(GRAPHS+'/*.fastg') + glob(GRAPHS+'/.gfa')
    reference = glob(GENOMES+'/reference.fasta')
    aligner = mp.Aligner(reference)
    genomes = list(aligner.seq_names)
    genome_lens = [len(aligner.seq(g)) for g in genomes]

    records = list()
    for graph in graphs:
        name = graph.replace('.fasta', '').replace('.gfa', '')
        asm = Assembly(assembler=name, assembly_file=graph)
        alignments, unaligned_contigs = align_nodes(asm, aligner)
        g_to_aln_idx = defaultdict(set)
        for i, aln in enumerate(alignments):
            g_to_aln_idx[aln.genome].add(i)
        contig_map=None
        cov_records = get_coverage_metrics(
            name, alignments, genomes, genome_lens, asm, g_to_aln_idx, unaligned_contigs, contig_map, disregard_unmapped=False
        )

        cov_records = cov_records.groupby('metric').value.sum()
        prc = cov_records.cov_tp / (cov_records.cov_tp + cov_records.cov_fp)
        rcl = cov_records.cov_tp / (cov_records.cov_tp + cov_records.cov_fn)
        fsc = (2*prc*rcl) / (prc+rcl)
        n_50 = n50(asm)
        nodes = num_nodes(asm)
        edges = num_edges(asm)
        param, asm = parse_name(graph) 
        records.append([asm, param, name, n_50, nodes, edges, 'cov_tp', cov_records.cov_tp])
        records.append([asm, param, name, n_50, nodes, edges, 'cov_fn', cov_records.cov_fn])
        records.append([asm, param, name, n_50, nodes, edges, 'cov_fp', cov_records.cov_fp])
        records.append([asm, param, name, n_50, nodes, edges, 'cov_precsion', prc])
        records.append([asm, param, name, n_50, nodes, edges, 'cov_recall', rcl])
        records.append([asm, param, name, n_50, nodes, edges, 'cov_F-score', fsc])


    output = pd.DataFrame(records, ['parameter', 'assembler', 'name', 'n50', 'num_nodes', 'num_edges', 'metric', 'value'])
    output.to_csv(os.path.join(OUTDIR, f'coverage_records.csv'), index=False)