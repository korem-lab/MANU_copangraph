#!/bin/bash 
#SBATCH --job-name=cami-gc
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/CAMISIMGraphQuality/log_%x-%j.log
python ./CAMISIMGraphQuality/run.py $ANALYSIS