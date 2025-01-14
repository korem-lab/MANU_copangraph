#!/bin/bash 
#SBATCH --job-name=momspi
#SBATCH --time=4:0:0
#SBATCH --mem=12G
#SBATCH --account=pmg
#SBATCH --cpus-per-task=4
#SBATCH --exclude=m013
#SBATCH --output=map-%x-%j.log
reads_dir=/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/HGF
contigs_dir=/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/EXT/megahit
final_dir=/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/MAP
# move to pmglocal
r1=$reads_dir/${SAMPLE}_R1.fastq.gz 
r2=$reads_dir/${SAMPLE}_R2.fastq.gz
out=/pmglocal/ic2465/MAP
mkdir $out
cp $r1 $r2 $out

#make index
bowtie2-build --threads 4 $contigs_dir/${SAMPLE}/final.contigs.fa $out/${SAMPLE}.idx 

# map reads
bowtie2 --threads 4 -x $out/${SAMPLE}.idx -1 $r1 -2 $r2 | samtools view -b -o $out/${SAMPLE}_mapping.bam
# count
samtools sort --threads 4 -o $out/${SAMPLE}_sorted_mapping.bam $out/${SAMPLE}_mapping.bam
samtools index --threads 4 -o $out/${SAMPLE}_sorted_mapping.bai $out/${SAMPLE}_sorted_mapping.bam
samtools idxstats $out/${SAMPLE}_sorted_mapping.bam > $out/${SAMPLE}_idxstats.csv
mv $out/${SAMPLE}*  $final_dir

