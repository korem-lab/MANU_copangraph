import os
import sys
import pandas as pd
import subprocess
from NodeBalancer import *
from glob import glob
import re

some_list = sorted(glob('/manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/analyses/homology_evaluation/scripts/input_yamls/*.yaml'))
some_list = [e for e in some_list if '10M' in e and '0.001' not in e]

#filt_list = list()
#for g in some_list:
#    sample = re.findall('sample_([0-9]+)_', g)[0]
#    asm = re.findall('_(megahit|metaspades)_', g)[0]
#    dt = re.findall('(0\.[0-9]+)\.yaml', g)[0]
#    if dt == '0.20':
#        dt = '0.2'
#    test_file = f'sample_{sample}_{asm}_10M_btwn_ani_avpv_{dt}.pkl'
#    test_file = f'sample_{sample}_{asm}_10M_alignment_objs_{dt}.pkl'
#    test_file =  os.path.join('/manitou/pmg/users/ic2465/homo_evaluation/results',test_file)
#    if not os.path.exists(test_file):
#        print(test_file)
#        filt_list.append(g)
#print(len(filt_list))
#    dump = [l.strip() for l in f]
#    filt_list = list()
#    for g in dump:
#        sample = re.findall('(sample_[0-9]+)_', g)[0]
#        asm = re.findall('_(megahit|metaspades)_', g)[0]
#        dt = re.findall('(0\.[0-9]+\.)pkl', g)[0]
#        selected = False
#        for e in some_list:
#            if sample in e and asm in e and dt in e:
#                selected = True
#                filt_list.append(e)
#        if not selected:
#            print(sample, asm, dt)
#

nb = NodeBalancer(limit=96)
num_cpu = 8
BOTTOM, TOP = 0, 50
for sample in some_list[BOTTOM:TOP]:
    node = nb.get_node(num_cpus=num_cpu)
    print('Fetched node %s!' % node)
    basename = os.path.basename(sample).replace('.yaml', '')
    fname = os.path.basename(sample)
    bash_file = f'launch_scripts/{basename}_homo_results.sh'

    bash_script = f"""#!/bin/bash 
#SBATCH --job-name={basename[-20:-10]}
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={num_cpu}
#SBATCH --account pmg
#SBATCH --output=/pmglocal/ic2465/_logs/%x-%j.{basename}.homo_analysis.log
#SBATCH --no-requeue
#SBATCH --nodelist={node}
#SBATCH --signal=B:TERM@60

echo '--------------------'
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST
echo '--------------------'

/usr/bin/time -v python /manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/analyses/homology_evaluation/scripts/homology_evaluation.py {sample}
"""

    with open(bash_file, "w") as file:
        file.write(bash_script)

#    subprocess.run(["sbatch", bash_file])
