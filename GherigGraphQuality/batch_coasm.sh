#!/bin/bash 
#SBATCH --job-name=coasm
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=32
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/coassemblies/log_%x-%j.log

r1=`cat $SAMPLES | xargs -n1 -I{} printf "/manitou/pmg/projects/korem_lab/Projects/Gherig/mmmbp/tmp/HGF2/%s_1.fastq.gz\n" {} | paste -sd,`
r2=`cat $SAMPLES | xargs -n1 -I{} printf "/manitou/pmg/projects/korem_lab/Projects/Gherig/mmmbp/tmp/HGF2/%s_2.fastq.gz\n" {} | paste -sd,`

echo $r1
echo $r2
megahit -t 32 -1 $r1 -2 $r2 -o $OUTDIR
