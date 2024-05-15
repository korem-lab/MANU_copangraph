#!/bin/bash 
#SBATCH --job-name=gherig-gc
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/log_%x-%j.log
python ./GherigGraphQuality/run.py $ANALYSIS
