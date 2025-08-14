#!/bin/bash 
#SBATCH --job-name=prokka
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/plt/prokka-%x-%j.log
mkdir -p /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/plt/TaskI/prokka/$SAMPLE
prokka --force --cpus 4 --outdir /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/plt/TaskI/prokka/ --prefix $SAMPLE /manitou/pmg/projects/korem_lab/Projects/ACU_PLT/mmmbp2/tmp_path2/EXT/megahit/$SAMPLE/final.contigs.fa

