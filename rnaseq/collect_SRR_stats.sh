#!/bin/bash 

TRIPLE_TABLE=$1

KK=`cut -f 1 $TRIPLE_TABLE`

for i in $KK
do
  PCT=`cat json/$i.run_info.json | grep -P "n_processed|n_pseudoaligned|p_pseudoaligned" | tr '[,\n]' '\t' | awk '{print $2"\t"$4"\t"$6}'`
  IDS=`grep -P "$i\t" $TRIPLE_TABLE`
  echo -e $IDS"\t"$PCT | tr ' ' '\t'
done
