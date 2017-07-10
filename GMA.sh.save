#!/bin/bash
#  Author: Marzie E. Rasekh
#  Laboratory for Biocomputing and Informatics, Boston University
#  contact me at marzie@bu.edu

## param
export PROC=8	                       # not using multi theaded for now
## sequencer errors
export SNP=0.01
export COVERAGE=1                     # total reads = (|sequence| - READ_LENGTH)*COVERAGE
export READ_LENGTH=101                 # %in% O(coverage*|sequence|*|read|) , at read length 100 coverage 1 for human chr1 (250M) makes 24 fastq files of size 2.5G each,
	                               # bam file sizes are neglectible (each ~200-250M)
export PADDING=100
export REFERENCE=/projectnb/vntrseek/share/GRCh38/GRCh38.fa
export BED_FILE=TRs.flanks.bed

## DIRECTORIES ##
export FASTA_DIR=fasta                 # to find $sequence.fasta file, notice the format
export FASTQ_DIR=fastq                 # to save/read the fastq files, suffix fq
export ALIGN_DIR=aln                   # to save/read alignment files
export FILE_BREAK=10000000             # setting this line break too low will add over-head, too high will make it slow
	                               # keep in mind that for each read in the fastq, there will be 4 lines
export DELIMINATOR="_"                 # this is set in the src/GMA_test.cpp code, %TODO fix this

USER=marzie
PROJECT=vntrseek


## initialize

if [ -d "$FASTA_DIR" ]; then
        rm $FASTA_DIR/*
else
    	mkdir $FASTA_DIR
fi

if [ -d "$FASTQ_DIR" ]; then
	rm $FASTQ_DIR/*
else
	mkdir $FASTQ_DIR
fi 

if [ -d "$ALIGN_DIR" ]; then
        rm $ALIGN_DIR/*
else
        mkdir $ALIGN_DIR 
fi 

## build
g++ src/GMA_SE.cpp -o bin/GMA_SE

module load samtools
module load bwa
module load bedtools
module load bedtools

awk -v DIR=$FASTA_DIR '{if(substr($1,1,1)==">") {F=substr($1,2);} print > DIR "/" F ".fasta"}' $REFERENCE

FILES=$FASTA_DIR/*.fa*
jobs=""
N=$(ls $FILES | wc -l)
for file in $FILES
do 
	i=$((i+1))
	qsub -u $USER -P $PROJECT -N GMA_seq$i -j y -pe omp $PROC -V GMA_per_fasta.sh $file
	jobs+=GMA_seq$i
	if [ "$i" -lt "$N" ]
	then
		jobs+=","
	fi
done

# while the jobs are not done wait
qsub -u $USER -P $PROJECT -hold_jid $jobs -l h_rt=100:00:00 -V GMA_TR_analysis.sh

