import os 
import sys
import pandas as pd 
import subprocess
from NodeBalancer import *
from glob import glob
from itertools import product

distances = [
    '0.005', '0.01', '0.02', '0.03',
    '0.04', '0.05', '0.06', '0.08', '0.1',
    '0.12', '0.14', '0.16', '0.18', '0.2', '0.25'
]

nb = NodeBalancer(limit=96)

num_cpu = 4
depth = [10]
combos = list(product(distances, range(10), ['megahit', 'metaspades'], depth))
print(len(combos))
BOTTOM, TOP = 250, 300
for dist, s, a, d in combos[BOTTOM:TOP]:
    sample = f'/manitou/pmg/users/ic2465/camisim_data/sample_{s}/{a}_sample_{d}M/sample_{s}_{a}_{d}M.pe_ext.fasta'
    node = nb.get_node(num_cpus=num_cpu, include=['m009', 'm013', 'm005',  'm004', 'm003', 'm006', 'm007'])
    print('Fetched node %s!'%node)
    basename = os.path.basename(sample).replace('.pe_ext.fasta', '')
    fname = os.path.basename(sample)
    bash_file = f'scripts/{basename}_{dist}_copan.sh'

    bash_script = f"""#!/bin/bash 
#SBATCH --job-name={s}-{a[:3]}-{d}-{dist}
#SBATCH --time=4:00:00
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --account pmg
#SBATCH --ntasks-per-node={num_cpu}
#SBATCH --output=/pmglocal/ic2465/_logs/%x-%j.sd={dist}.log
#SBATCH --no-requeue
#SBATCH --nodelist={node}

echo '--------------------'
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST
echo '--------------------'

/usr/bin/time -v python /manitou/pmg/users/ic2465/pg_collapse_algorithm/src/main.py -c {sample} -g {basename}_sd{dist} -o /manitou/pmg/users/ic2465/copan_manuscript/homo_evaluation/graph_out_nosm -d {dist} --asymmetrical --threads {num_cpu}
    """
    with open(bash_file, "w") as file:
        file.write(bash_script)

    #subprocess.run(["sbatch",bash_file])
