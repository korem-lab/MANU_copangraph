#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=64
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/PanMetagenomeViz/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/PanMetagenomeViz
gtdbtk=/pmglocal/ic2465/gtdbtk

mkdir ${gtdbtk}
mkdir ${gtdbtk}/genomes
cp ${data}/*bin*.fa ${gtdbtk}/genomes
gtdbtk classify_wf --genome_dir ${gtdbtk}/genomes --out_dir ${gtdbtk}/classify -x fa --cpus 64
cp -r ${gtdbtk}/classify ${data}

