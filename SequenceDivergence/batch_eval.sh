#!/bin/bash 
#SBATCH --job-name=seq_div
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/SequenceDivergence/log_%x-%j.log
python ./SequenceDivergence/run_sequence_divergence.py $ANALYSIS
