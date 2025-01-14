#!/bin/bash 
#SBATCH --job-name=ms
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/crc/extract-%x-%j.log
python extract_subgraph_from_gfa.py CRC_CN $TOPFS 10 ../data/Predictions/crc/ ../data/Predictions/crc/CRC_CN_sample_names.list
