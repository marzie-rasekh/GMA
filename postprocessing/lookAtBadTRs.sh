#!/bin/bash

### silly code to get sequences and flanks of bad TRs ###

# give me those bad TRs,format: ID chrom start end then whatever
TR_LIST_FILE="$1"
# give me the flank size of interest
FLANK_LENGTH="$2"
OFFSET="$3"

REFERENCE=/projectnb/vntrseek/share/GRCh38/GRCh38.fa

sort -n -k 5 $TR_LIST_FILE > tmp
mv tmp $TR_LIST_FILE

# this will produce output as ID chrom start end whatever flank_left flank_right TR

echo "reading " $TR_LIST_FILE 
rm test
cat $TR_LIST_FILE | while read line
do
   echo $line | awk -v R=$REFERENCE -v F=$FLANK_LENGTH -v O=$OFFSET '{print "samtools faidx " R " "  $2 ":" $3-F-1+O "-" $3-1+O}' > tmp.sh
   LEFT=$(sh tmp.sh | awk '(NR>1) {print}' | tr -d "\n")
   echo $line | awk -v R=$REFERENCE -v F=$FLANK_LENGTH -v O=$OFFSET '{print "samtools faidx " R " "  $2 ":" $4+1-O "-" $4+F+1-O}' > tmp.sh
   RIGHT=$(sh tmp.sh | awk '(NR>1) {print}' | tr -d "\n")
   echo $line | awk -v R=$REFERENCE -v O=$OFFSET '{print "samtools faidx " R " "  $2 ":" $3+O "-" $4-O}' > tmp.sh
   ARRAY=$(sh tmp.sh | awk '(NR>1){print}' | tr -d "\n")
   echo $line$'\t'$LEFT $RIGHT $ARRAY >> test
done
rm tmp.sh

awk '{print $1 "\t" $2 ":" $3 "-" $4 "\t" ($12 ? "I" : "S") "\t" ($13 ? "-" : "M") "\t" ($4-$3+1) "\t"  $14 " " $15 " " $16}' test | sort -f -k 5,5 -k 6,6 | uniq -f 5 -i -D > unmarked.indis.BADHOMBRE.txt
rm test
