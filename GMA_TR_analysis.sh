#!/bin/bash
#  Author: Marzie E. Rasekh
#  Laboratory for Biocomputing and Informatics, Boston University
#  contact me at marzie@bu.edu

cat $ALIGN_DIR/*.extra.bed > $ALIGN_DIR/extra.bed
cat $ALIGN_DIR/*.missing.bed > $ALIGN_DIR/missing.bed

bedtools coverage -counts -f 1.00 -a $BED_FILE -b $ALIGN_DIR/missing.bed > $ALIGN_DIR/missing.coverage
bedtools coverage -counts -f 1.00 -a $BED_FILE -b $ALIGN_DIR/extra.bed   > $ALIGN_DIR/extra.coverage
paste $ALIGN_DIR/missing.coverage $ALIGN_DIR/extra.coverage | awk -v R=$READ_LENGTH '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $3-$2 "\t" $5 "\t" $10 "\t" (R-($3-$2)+1-$5)/(R-($3-$2)+1+$10)}' > $ALIGN_DIR/calculated.gma

rm $ALIGN_DIR/missing.coverage
rm $ALIGN_DIR/extra.coverage
