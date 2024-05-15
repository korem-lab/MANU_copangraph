
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
    complexity = pd.DataFrame(columns=['assembler', 'n_samples', 'rep', 'ms', 'metric', 'value'])
    # populate complexity table
    for gfa in gfas:
        n_samples, rep, ms = re.findall('coasm_([0-9]+)_rep_([0-9])_ms_([0-9]+).gfa', gfa)[0]
        asm = Assembly('copangraph', gfa)
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'nodes', num_nodes(asm)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'edges', num_edges(asm)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'N50', compute_nX(asm.contigs, 50)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'N90', compute_nX(asm.contigs, 90)]
    for n_samples, rep in product(N_SAMPLES, REPS):
        fastg = glob.glob(os.path.join(DATA, f'coasm_{n_samples}_rep_{rep}/*.fastg'))[0]
        asm = Assembly('megahit', fastg)
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'nodes', num_nodes(asm)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'edges', num_edges(asm)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'N50', compute_nX(asm.contigs, 50)]
        complexity.loc[len(complexity), :] = [asm.assembler, n_samples, rep, ms, 'N90', compute_nX(asm.contigs, 90)]
    complexity.to_csv(os.path.join(DATA, 'complexity_results.csv'))