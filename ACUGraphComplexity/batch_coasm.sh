#!/bin/bash 
#SBATCH --job-name=ms
#SBATCH --time=8:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/log_%x-%j.log
# assemble
r1=`cat /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.txt | xargs -I {} printf "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/fastq/%s_1.fastq.gz," {}`
r2=`cat /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.txt | xargs -I {} printf "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/fastq/%s_2.fastq.gz," {}`
/usr/bin/time -v megahit -1 $r1 -2 $r2 -t 32 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/coassembly_${REP}



