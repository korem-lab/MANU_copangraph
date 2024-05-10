#!/bin/bash 
#SBATCH --job-name=ms
#SBATCH --time=8:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=64
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity
write=/pmglocal/ic2465


# pool
cat $SAMPLES | while read s; do
    s=`echo $s | rev | cut -b 13`
	echo ${s} ${COASSEMBLY}
	zcat $data/${s}_R1.fastq.gz >> ${write}/${COASSEMBLY}_pooled_R1.fastq
	zcat $data/${s}_R2.fastq.gz >> ${write}/${COASSEMBLY}_pooled_R2.fastq
done

# compress
pigz -p 64 ${write}/${COASSEMBLY}_pooled_R1.fastq
pigz -p 64 ${write}/${COASSEMBLY}_pooled_R2.fastq

# assemble
/usr/bin/time -v metaspades -1 ${write}/${COASSEMBLY}_pooled_R1.fastq.gz -2 ${write}/${COASSEMBLY}_pooled_R2.fastq.gz -t 64 -o ${write}/${COASSEMBLY}

# copy
mv ${write}/${COASSEMBLY} ${data}
mv ${write}/${COASSEMBLY}_pooled_R[12].fastq.gz ${data}


