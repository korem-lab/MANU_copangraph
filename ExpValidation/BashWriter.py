import os 
import pandas as pd 
import subprocess
from NodeBalancer import *
from pathlib import Path
from os.path import join
import uuid
import sys
import shutil
from functools import partial
import re

class BashWriter():
    def __init__(self,
        mem=12, # mem (in Gb) requested per job
        hours=24, # number of hours requested per job
        cpus=2, # number of CPUs requested per job
        out_dir=None, # directory where logs/ and scripts/ files are created (default to cwd)
        include=None, # list of included nodes (keep None unless you want to target some specific nodes)
        exclude=None, # list of excluded nodes (keep None unless you want to exclude some specific nodes)
        limit=999, # limit number of CPUs on each node
        storage=None # minimum % of storage needed on /pmglocal/ node
    ):
        self.mem = int(mem)
        self.hours = int(hours)
        self.cpus = int(cpus)
        self.storage = int(storage) if storage is not None else None
        self.batch = []

        out_dir = os.getcwd() if out_dir is None else out_dir
        self.log_dir = os.path.join(out_dir,'logs')
        self.script_dir = os.path.join(out_dir,'scripts')
        
        self.user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
        self.nb = NodeBalancer(limit=limit,include=include,exclude=exclude)
       
        if os.path.isdir(self.script_dir):
            shutil.rmtree(self.script_dir)

        print('Log directory- %s'%self.log_dir)
        print('Script directory- %s'%self.script_dir)

        Path(self.log_dir).mkdir(parents=True, exist_ok=True)
        Path(self.script_dir).mkdir(parents=True, exist_ok=True)
  
    def add(self,command):
        self.batch.append(command)

    def run_batch(self,cmds_per_job=5,name='bw',dry=False,remove=True,check_if_running=False,stagger_seconds=0,out_dir=None,sleep=0,parallel=False,remove_copied_files=True):
        groups = list(_split(self.batch,cmds_per_job))
        print('Running %s jobs with %s commands each...'%(len(groups),cmds_per_job))
        for group in groups:
            command = '\n'.join(group) if not parallel else ' &\n'.join(group) + ' &\nwait'
            self.run(command,name,dry,remove,check_if_running,stagger_seconds,sleep,out_dir,remove_copied_files)
        self.batch = []

    def run(self,
        command, # command you want to run
        name='bw', # job name will be this variable followed by random uuid4
        dry=False, # set this to True if you just want to create the bash script but not actually run
        remove=True, # set this to True if you want to delete each bash script after submitting with `sbatch`
        check_if_running=False, # set this to True if you want to check to see if a job is already running (based off of `name` parameter)
        stagger_seconds=0, # Sleep random number of seconds before beginning job (help stagger load on IO heavy jobs)
        sleep=0, # Sleep set number of seconds before beginning jobs
        out_dir=None, # Directory you want to move output files into from /pmglocal/ directory
        remove_copied_files=True # Remove any files you copied over into /pmglocal/
    ):
        if check_if_running:
            curr_jobs = [x[:-6] for x in get_jobs()]
            if name in curr_jobs:
                print('%s already running! Skipping...'%name)
                return
        job_name = f"{name}_{str(uuid.uuid4())[:5]}"
        bash_file = os.path.join(self.script_dir,job_name + '.sh')
        
        command = '\n# --START --\n' + '\n'.join([c.strip() for c in command.split('\n') if len(c.strip())>0]) + '\n# --END --\n\n'

        sleep_cmd = ""
        if stagger_seconds > 0:
            sleep_cmd = f"\nMINWAIT=10\nMAXWAIT={stagger_seconds}\nRANDOM_WAIT=$((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))\necho Sleeping $RANDOM_WAIT seconds...\nsleep $RANDOM_WAIT"
        elif sleep > 0:
            sleep_cmd = f"echo Sleeping {sleep} seconds...\nsleep {sleep}"
        pmglocal_run_dir = f"/pmglocal/{self.user}/bashwriter/{job_name}"
        files_to_copy = re.findall(r'\$\$(.*?)\$\$', command)
        
        if ('$CUID$' in command):
            command = command.replace('$CUID$',self.user)

        if ('$PMGLOCAL$' in command or len(files_to_copy) > 0):
            for file in files_to_copy:
                command = f"cp -r {file} $PMGLOCAL$/ && echo Copied {file} to run directory...\n" + command
                command = command.replace(f"$${file}$$",f"$PMGLOCAL$/{os.path.basename(file)}")
                if remove_copied_files:
                    command = command + f"rm -f $PMGLOCAL$/{os.path.basename(file)}\n"
            command = f"mkdir $PMGLOCAL$/tmp/ && echo Created temporary directory $PMGLOCAL$/tmp/ && export TMPDIR=$PMGLOCAL$/tmp\n" + command
            command = f"mkdir -p $PMGLOCAL$ && echo Created /pmglocal/ run directory $PMGLOCAL$...\ncd $PMGLOCAL$\n" + command
            command = command.replace('$PMGLOCAL$',pmglocal_run_dir)
            if out_dir is not None:
                command = command + f"\necho '\nMoving output to {out_dir}...'\nmkdir -p {out_dir}\nmv {pmglocal_run_dir}/* {out_dir}/ && echo Done!"
            else:
                print('!!WARNING!! `out_dir` set to None! Output files will not be saved...')
            command = command + f"\nrm -rf {pmglocal_run_dir}/"
        
        node = self.nb.get_node(num_cpus=self.cpus,storage_space=self.storage)
        
        scaffold = f"""#!/bin/bash 
#SBATCH --job-name={job_name}
#SBATCH --time={self.hours}:00:00
#SBATCH --mem={self.mem}G
#SBATCH --ntasks-per-node={self.cpus}
#SBATCH --account pmg
#SBATCH --output={self.log_dir}/%j-%x.log
#SBATCH --nodelist={node}
#SBATCH --no-requeue

echo '--------------------'
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST
echo '--------------------'
{sleep_cmd}
{command}
    """
        with open(bash_file, "w") as file:
            file.write(scaffold)
        if not dry:
            subprocess.run(["sbatch",bash_file]) 
        if remove:
            os.remove(bash_file)

def _split(list_a, chunk_size):
    for i in range(0, len(list_a), chunk_size):
        yield list_a[i:i + chunk_size]
