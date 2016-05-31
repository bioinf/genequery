#!/bin/bash

# Converts GPLXXX_family.soft_.gz file into 3-column annotation table: 1) probe 2) symbol 3) entrez
#
# Writes GPL###.3col file to current directory as a result.

if [ "$#" -ne 4 ]; then
    echo "USAGE: bash soft_to_3col.sh <path to GPLXXX_family.soft.gz> <probe id column> <symbol id column> <entrez id column>"
    exit 1
fi

PATH_TO_GPL=$1
GPL=${PATH_TO_GPL%%_*}
GPL=${GPL##*/}

P=$2
S=$3
E=$4

## all absent fields should be replaced with NONE string

zcat $PATH_TO_GPL | sed -n '/!platform_table_begin/,/!platform_table_end/ p' | sed '1d;$d' | tail -n +2 | awk -F "\t" -v a=$P -v b=$S -v c=$E '{print $a "\t" $b "\t" $c }' | sed 's/\t\t/\tNONE\t/g' | sed 's/\t$/\tNONE/g' > ${GPL}.3col
