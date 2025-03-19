#!/bin/bash 
#SBATCH --job-name=gherig-mixing
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=m005
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/repmix-%x-%j.log
python /burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/gherig_run.py $ANI
