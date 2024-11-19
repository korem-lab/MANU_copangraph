#!/bin/bash 
#SBATCH --job-name=cnx_cov
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=%x-%j.log
cd /burg/pmg/users/ic2465/Projects/MANU_copangraph
python CAMISIMGraphQuality/run.py ./data/CAMISIMCoassembly/coasm_${SAMPLE}_analysis.yaml
