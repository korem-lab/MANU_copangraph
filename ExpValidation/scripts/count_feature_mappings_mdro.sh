#!/bin/bash 
#SBATCH --job-name=cntfeat
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=m013
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ExpValidation/feature_counts/mdro/fcounts-%x-%j.log
FEATURE=$1
fasta=/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/Predictions/mdro/scratch/copangraphs/mdro_pos_02_mo1000_ms75.fasta
ncol=/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/Predictions/mdro/scratch/copangraphs/mdro_pos_02_mo1000_ms75.ncolor.csv
ecol=/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/Predictions/mdro/scratch/copangraphs/mdro_pos_02_mo1000_ms75.ecolor.csv
outdir=/pmglocal/ic2465
read_names=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ExpValidation/feature_counts/mdro/dummy.list
samples=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ExpValidation/feature_counts/mdro/mdro_pos.list
persistence_table=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ExpValidation/feature_counts/mdro/persistence_table.csv

# make feature speicifc directory and clean up prior analysis directories
outdir=$outdir/$FEATURE
#rm -rf $outdir
#mkdir $outdir

# extract feature fasta from graph fasta
# a binary vector of which samples contribute to the feature
# and a vector of the lengths of each sequence
#python extract_feature_fasta.py $FEATURE $fasta $ncol $ecol $samples $outdir

# index the feature
#bowtie2-build $outdir/${FEATURE}.fasta $outdir/feature

# map reads to the features and count the mappings
cat $read_names | while read s; do
	reads_1=${s}_1.fastq.gz
	reads_2=${s}_2.fastq.gz
	echo bowtie
	bowtie2 --local --no-unal -a --threads 16 -x $outdir/feature -U $reads_1 -S ${outdir}/mapping_1.sam
	bowtie2 --local --no-unal -a --threads 16 -x $outdir/feature -U $reads_2 -S ${outdir}/mapping_2.sam
	echo samtools
	samtools sort -@ 16 -n ${outdir}/mapping_1.sam -o ${outdir}/mapping_1.sort.sam
	samtools sort -@ 16 -n ${outdir}/mapping_2.sam -o ${outdir}/mapping_2.sort.sam
	echo rm ${outdir}/mapping.sam
	rm ${outdir}/mapping_1.sam
	rm ${outdir}/mapping_2.sam
	echo counting

	# get num reads 
	num_reads=$(zcat "$reads_1" | wc -l)
	num_reads=$((num_reads / 2))  # divide by 4, and then times by two because its paired (so just divide by 2)
	echo $num_reads
	# produces a counter pickle file that counts how many reads map to each "reference sequence"
	# multimapping across the refs is allowed: a read can map to ref seq a, b, and both of their
	# counters will be incremented, but more than one mapping of a read to a will only increment the counter once. 
	# also outputs a total library size to normalize for this
	python tally_mappings.py ${outdir}/mapping_1.sort.sam ${outdir}/mapping_2.sort.sam ${s} $num_reads ${outdir}
done
## compute feature count vector
python construct_feature_count_table.py $FEATURE $outdir $persistence_table ${outdir}/${FEATURE}.occ.csv ${outdir}/${FEATURE}.seqlens.csv
mv $outdir /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ExpValidation/feature_counts/mdro
