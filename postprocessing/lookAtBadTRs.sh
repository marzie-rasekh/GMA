#!/bin/bash

### silly code to get sequences and flanks of bad TRs ###

# give me those bad TRs,format: ID chrom start end then whatever
TR_LIST_FILE="$1"
# give me the flank size of interest
FLANK_LENGTH="$2"
OFFSET="$3"

REFERENCE=/projectnb/vntrseek/share/GRCh38/GRCh38.fa

# this will produce output as ID chrom start end ... flank_left flank_right array

#echo "reading " $TR_LIST_FILE 
> test
line=0
cat $TR_LIST_FILE | while read id chrom start end size missing extra gma_score missing_normal extra_normal coverage indis multi
do
   if [ "$line" -eq "0" ]; then
      line=1
   else
      LEFT=$(samtools faidx $REFERENCE $chrom:$((start - FLANK_LENGTH + OFFSET))-$((start + OFFSET)) | awk '(NR>1) {print}' | tr -d "\n")
      RIGHT=$(samtools faidx $REFERENCE $chrom:$((end - OFFSET))-$((end + FLANK_LENGTH - OFFSET)) | awk '(NR>1) {print}' | tr -d "\n")
      ARRAY=$(samtools faidx $REFERENCE $chrom:$((start + OFFSET))-$((end - OFFSET)) | awk '(NR>1) {print}' | tr -d "\n")
      echo $id$'\t'$chrom:$start-$end$'\t'$size$'\t'$gma_score$'\t'$indis$'\t'$multi$'\t'$LEFT $RIGHT $ARRAY >> test
   fi
done

awk '{print $1 "\t" $2 "\t" $3 "\t" substr(100*$4,1,index(100*$4,".")-1) "." substr(100*$4,index(100*$4,".")+1,3) "%\t" ($5 ? "I" : "S") "\t" ($6 ? "-" : "M") "\t"  $7 " " $8 " " $9 }' test | sort -f -k 3,3 -k 7,7 | uniq -f 6 -i -D > ${TR_LIST_FILE%.*}.repeats 
rm test
