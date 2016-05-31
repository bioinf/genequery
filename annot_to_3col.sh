#!/bin/bash

## Converts GPLXXX.annot.gz file into 3-column annotation table: 1) probe 2) symbol 3) entrez

ANNOT=$1

GPL=${ANNOT%%.annot.gz}

if [[ $ANNOT == "" ]] 
then
  echo "ERROR: No input annotation file."
  exit 1
fi

if [ ! -s $ANNOT ]; then
   echo "$ANNOT has zero size"
   exit 1
fi

## all absent fields should be replaced with NONE string

zcat $ANNOT | sed -n '/!platform_table_begin/,/!platform_table_end/ p' | sed '1d;$d' | awk -F "\t" '{print $1 "\t" $3 "\t" $4}' | sed 's/\t\t/\tNONE\t/g' | sed 's/\t$/\tNONE/g' > $GPL.3col
