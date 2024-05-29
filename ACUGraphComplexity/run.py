
import os
import re
import pandas as pd
import glob
from itertools import product
from Utils.AssemblyParser import Assembly
from Utils.evaluate_assembly import num_nodes, num_edges, compute_nX

DATA = './data/ACUGraphComplexity'
N_SAMPLES = [1, 5, 10, 20, 50]
REPS = [0, 1, 2]

if __name__ == '__main__':


    gfas = glob.glob(os.path.join(DATA, 'coasm*.gfa'))
    complexity = pd.DataFrame(columns=['assembler', 'n_samples', 'rep', 'ms', 'dt', 'metric', 'value'])
    gfas_dt = [e for e in gfas if 'dt' in e]
    gfas = [e for e in gfas if 'dt' not in e]
    # populate complexity table
    #for gfa in gfas:
    #    n_samples, rep, ms = re.findall('coasm_([0-9]+)_rep_([0-9])_ms_([0-9]+).gfa', gfa)[0]
    #    asm = Assembly('copangraph', gfa)
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 0.02, 'nodes', num_nodes(asm)]
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 0.02, 'edges', num_edges(asm)]
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 0.02, 'N50', compute_nX(asm.contigs, 50)]
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 0.02, 'N90', compute_nX(asm.contigs, 90)]
    #for n_samples, rep in product(N_SAMPLES, REPS):
    #    fastg = glob.glob(os.path.join(DATA, f'coasm_{n_samples}_rep_{rep}/*.fastg'))[0]
    #    asm = Assembly('megahit', fastg)
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, pd.NA, 'nodes', num_nodes(asm)]
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, pd.NA, 'edges', num_edges(asm)]
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, pd.NA, 'N50', compute_nX(asm.contigs, 50)]
    #    complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, pd.NA, 'N90', compute_nX(asm.contigs, 90)]
    #complexity.to_csv(os.path.join(DATA, 'complexity_results.csv'))
    complexity = pd.DataFrame(columns=['assembler', 'n_samples', 'rep', 'ms', 'dt', 'metric', 'value'])
    for gfa in gfas_dt:
        n_samples, rep, ms, dt = re.findall('coasm_([0-9]+)_rep_([0-9])_ms_([0-9]+)_dt_(0.[0-9]+).gfa', gfa)[0]
        asm = Assembly('copangraph', gfa)
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, dt, 'nodes', num_nodes(asm)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, dt, 'edges', num_edges(asm)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, dt, 'N50', compute_nX(asm.contigs, 50)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, dt, 'N90', compute_nX(asm.contigs, 90)]
    complexity.to_csv(os.path.join(DATA, 'complexity_results_dt.csv'))