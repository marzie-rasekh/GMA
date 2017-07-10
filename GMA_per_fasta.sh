#!/bin/bash

file="$1"
sequence=$(basename "$file" | cut -d. -f1)			# the sequence to process, only the name, will find FASTA_DIR/{$sequence}.fa*

echo $sequence

## make sample reads
bin/GMA_SE $FASTA_DIR/$sequence.fa* -SNP_RATE $SNP -COVERAGE $COVERAGE -READ_LENGTH $READ_LENGTH -FILE_BREAK $FILE_BREAK


if [ -e $ALIGN_DIR/$sequence.gma ]
then
    rm $ALIGN_DIR/$sequence.gma
fi

# map the reads and extract the nonmappers
FILES=$FASTQ_DIR/$sequence$DELIMINATOR*.fq
for fq_file in $FILES
do
    ## map the reads
    echo "mapping and processing " $fq_file 
    ## extract the nonmapping reads as : 
    ## original_sequence original_left_position mapped_sequence mapped_left_position
    ## given the format of the read IDs are : se_[sequence]_[leftposition]
    ## only works for sequences of max length 999999999, if you need larger please change the padding in the src/GMA_test.cpp code (currently 10) 
    ## extract fields and check non-mapping by padding
    bwa mem -t $PROC $REFERENCE $fq_file | awk '{if(substr($1,1,1)!="@") print substr($1,4,index($1,"_0")-4) "\t" (0+(substr($1, index($1,"_0") + 1, index($1,"(") - index($1,"_0") - 1) )) "\t" $3 "\t" $4 "\t" ($3=="*"?"*":(int(($2%32)/16)==1?"+":"-")) }' |  awk -v PADD=$PADDING '{if($1!=$3 || !((($4-$2) < PADD) && (($4-$2) > (-1*PADD)))) print }' >> $ALIGN_DIR/$sequence.gma
    ## clean up to save space
    rm $fq_file
done

# make the missing and extra files
awk -v R=$READ_LENGTH '{C=0; if (N==$2) {N=$2; C++;} else {N=$2; print $1 "\t" $2 "\t" ($2+R) "\t" $5 "\t" C; C=0;}}' $ALIGN_DIR/$sequence.gma | bedtools sort -i - > $ALIGN_DIR/$sequence.missing.bed
awk -v R=$READ_LENGTH '{if($3!="*") {print $3 "\t" $4 "\t" ($4+R) "\t" $5 "\t"} }' $ALIGN_DIR/$sequence.gma | bedtools sort -i - > $ALIGN_DIR/$sequence.extra.bed







