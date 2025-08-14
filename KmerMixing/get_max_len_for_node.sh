#!/bin/bash 
#SBATCH --job-name=maxnode
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=m004
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/maxnode-%x-%j.log
python KmerMixing/get_max_len_for_node.py $GRAPH 
