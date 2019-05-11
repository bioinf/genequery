#!/bin/bash 

GSM_LIST=$1
GSE_SEL=$3
GSE_META=$3
SPECIES=$4

if [[ $# != 4 ]]
then
  echo "USAGE: ./select_best_gse2.sh Human_GSM.list Human_GSE_selected_v2.tsv Human_GSE_metadata_v2.tsv \"Homo sapiens\""
  exit 1
fi 

KK=`cat $GSM_LIST`

for i in $KK
do
  PP=`grep -P "\t$i\t" $GSE_SEL | cut -f 1`
  a=( $PP ) 
  if [[ ${a[1]} != "" ]]
  then
    BEST_PCT=0
    BEST_GSE=""
    SIZE_GSE=""

    >&2 echo "GSM entry $i has MORE than one GSE: $PP"
    for j in $PP 
    do
      PCT=`grep -P "$j\t" $GSE_META | awk -v v="$SPECIES" -F "\t" '{if ($3=="cDNA" && $4==v && $5=="RNA-Seq") {sum++}} END {printf "%d\n",100*sum/NR}'`
      >&2 echo "GSE $j contains $PCT% RNA-Seq experiments"
      if (($PCT > $BEST_PCT)) 
      then 
        BEST_PCT=$PCT
        BEST_GSE=$j
      fi
    done
    >&2 echo "Final decision for GSM $i: best GSE is $BEST_GSE"
    echo -e $i"\t"$BEST_GSE 
  else 
    # >&2 echo "GSM entry $i matches one GSE: $PP"
    echo -e $i"\t"$PP
  fi
done
