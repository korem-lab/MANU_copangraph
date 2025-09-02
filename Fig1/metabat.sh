#!/bin/bash 
#SBATCH --job-name=bt2
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/PanMetagenomeViz/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/PanMetagenomeViz
write=/pmglocal/ic2465/pmgviz_${SAMPLE}

rm -rf ${write}
mkdir ${write}
cp ${data}/${SAMPLE}.final.contigs.fa ${write}
cp ${data}/${SAMPLE}.bam ${write}
jgi_summarize_bam_contig_depths --outputDepth ${write}/${SAMPLE}.depth.txt ${write}/${SAMPLE}.bam
metabat2 -i ${write}/${SAMPLE}.final.contigs.fa -a ${write}/${SAMPLE}.depth.txt -o ${write}/${SAMPLE}.bin
cp ${write}/${SAMPLE}.bin*.fa ${data}
rm -rf ${write}



