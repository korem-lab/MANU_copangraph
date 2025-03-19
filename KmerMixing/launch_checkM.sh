#!/bin/bash 
#SBATCH --job-name=checkM
#SBATCH --time=4:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --exclude=m005
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/checkM-%x-%j.log
checkm lineage_wf -x fa $DIR ${DIR}/checkm_out -t 16 --file output.txt
