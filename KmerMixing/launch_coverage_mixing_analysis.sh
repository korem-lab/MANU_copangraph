#!/bin/bash 
#SBATCH --job-name=gherig-mixing
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/repmix-%x-%j.log
/usr/bin/time -v python KmerMixing/mixing_calcs.py $GRAPH $WINDOW
