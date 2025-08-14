#!/bin/bash 
#SBATCH --job-name=grg-ss
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --account pmg
#SBATCH --exclude=m005
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/grgss_%x-%j.log

TOTAL=80
tmp=/pmglocal/ic2465/${SAMPLE}_${REQUIRED}
rm -rf $tmp
mkdir $tmp

cp /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_1.fastq.gz ${tmp}/${SAMPLE}_1.fastq.gz
cp /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_2.fastq.gz ${tmp}/${SAMPLE}_2.fastq.gz
gunzip ${tmp}/*.gz

OUT1=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_${REQUIRED}M_1.fastq
OUT2=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_${REQUIRED}M_2.fastq

TMPOUT1=${tmp}/${SAMPLE}_${REQUIRED}M_1.fastq
TMPOUT2=${tmp}/${SAMPLE}_${REQUIRED}M_2.fastq
R1=${tmp}/${SAMPLE}_1.fastq
R2=${tmp}/${SAMPLE}_2.fastq
python -c "import parse_seq as ps; ps.subsample_reads('$R1', '$R2', '$TMPOUT1', '$TMPOUT2', $TOTAL, $REQUIRED, $SAMPLE)"

pigz -p 16 $TMPOUT1
pigz -p 16 $TMPOUT2
mv ${TMPOUT1}.gz /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_${REQUIRED}M_1.fastq.gz
mv ${TMPOUT2}.gz /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_${REQUIRED}M_2.fastq.gz
mv ${TMPOUT1}.seed /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_${REQUIRED}M_1.fastq.seed
mv ${TMPOUT2}.seed /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/${SAMPLE}_out/data/reads/${SAMPLE}_${REQUIRED}M_2.fastq.seed

rm -rf $tmp
