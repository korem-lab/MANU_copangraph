#!/bin/bash 
#SBATCH --job-name=ms
#SBATCH --time=8:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=64
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/
write=/pmglocal/ic2465


# clear
rm ${write}/${COASSEMBLY}_pooled_R1.fastq
rm ${write}/${COASSEMBLY}_pooled_R2.fastq
rm ${write}/${COASSEMBLY}_pooled_R1.fastq.gz
rm ${write}/${COASSEMBLY}_pooled_R2.fastq.gz

# pool
cat $SAMPLES | while read s; do
    s=`echo $s | rev | cut -b 14- |rev`
	echo ${s} ${COASSEMBLY}
	zcat ${s}_R1.fastq.gz >> ${write}/${COASSEMBLY}_pooled_R1.fastq
	zcat ${s}_R2.fastq.gz >> ${write}/${COASSEMBLY}_pooled_R2.fastq
done

# compress
pigz -p 64 ${write}/${COASSEMBLY}_pooled_R1.fastq
pigz -p 64 ${write}/${COASSEMBLY}_pooled_R2.fastq

# assemble
/usr/bin/time -v megahit -1 ${write}/${COASSEMBLY}_pooled_R1.fastq.gz -2 ${write}/${COASSEMBLY}_pooled_R2.fastq.gz -t 64 -o ${write}/${COASSEMBLY}

# copy
mv ${write}/${COASSEMBLY} ${data}
mv ${write}/${COASSEMBLY}_pooled_R[12].fastq.gz ${data}


