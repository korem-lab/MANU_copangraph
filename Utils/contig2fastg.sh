#!/bin/bash 
#SBATCH --job-name=grg-ss
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/Utils/contig-2fastg-%x-%j.log

megahit_toolkit contig2fastg $K $CONTIGS > $FASTG
