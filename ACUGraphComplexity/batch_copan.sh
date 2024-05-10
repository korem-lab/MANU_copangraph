#!/bin/bash 
#SBATCH --job-name=cpg
#SBATCH --time=1:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity
write=/pmglocal/ic2465
COPANGRAPH=/burg/pmg/users/ic2465/copangraph/bin/release/copangraph

echo ./usr/bin/time -v $COPANGRAPH $INI 
