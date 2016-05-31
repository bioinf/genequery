#!/bin/bash 

SPECIES=$1
## this is the script intended to make sure there are no redundant GSEs - such as super-series or simply redundant GEO submissions. 

for i in *_eigengenes.txt
do 
  GSE=${i%%_eigengenes.txt}
  head -n 1 $i | sed "s/\t/\n/g" | perl -sne 'm/(GSM\d+)/; print "$gse\t$1\n"' -- -gse=$GSE
done  > $SPECIES.gse_gsm.table 

awk '{print $2}' $SPECIES.gse_gsm.table | sort | uniq -c | sort -k1,1n > $SPECIES.gsm_freq.table
echo "Done generating GSM frequency table."

N2=`grep -c -P "2 GSM" $SPECIES.gsm_freq.table`
NN=`awk '{if ($1>2) print $2}' $SPECIES.gsm_freq.table | wc -l`

echo "Found $N2 GSMs that occur twice and $NN GSMs that occur more than twice."

## get GSMs that are encountered in exactly 2 GSEs 
KK1=`awk '{if ($1==2) print $2}' $SPECIES.gsm_freq.table`

for i in $KK1
do
  GSE1=`grep $i $SPECIES.gse_gsm.table | awk '{print $1}' | head -n 1`
  GSE2=`grep $i $SPECIES.gse_gsm.table | awk '{print $1}' | tail -n 1`
  N1=`grep -c -P "$GSE1" $SPECIES.gse_gsm.table`
  N2=`grep -c -P "$GSE2" $SPECIES.gse_gsm.table`
  N12=`grep -P "$GSE1|$GSE2" $SPECIES.gse_gsm.table | awk '{print $2}' | sort | uniq | wc -l`
  if [[ $((N1+N2)) > $N12 ]]
  then 
    echo -e "$GSE1\t$GSE2\t$((N1+N2-N12))\t$N1\t$N2"
  fi
done | sort | uniq > $SPECIES.gse_x2_overlap.table

echo "Done making a list of overlapping GSEs."
awk '{if ($4==$5) {sum1++} else {sum2++}} END {printf "Same size GSEs: %d, different size GSEs: %d, overall: %d.\n",sum1,sum2,sum1+sum2}' $SPECIES.gse_x2_overlap.table

## get GSMs that are encountered in >2 GSEs 
KK2=`awk '{if ($1>2) print $2}' $SPECIES.gsm_freq.table`

for i in $KK2 
do 
  grep -P "\t$i$" $SPECIES.gse_gsm.table | awk '{print $1}' 
done | sort | uniq > $SPECIES.gse_xN.list

for i in `cat $SPECIES.gse_xN.list`
do  
  for j in `cat $SPECIES.gse_xN.list`
  do
    if [[ $i < $j ]]
    then
      GSE1=$i
      GSE2=$j
      N1=`grep -c -P "$GSE1" $SPECIES.gse_gsm.table`
      N2=`grep -c -P "$GSE2" $SPECIES.gse_gsm.table`
      N12=`grep -P "$GSE1|$GSE2" $SPECIES.gse_gsm.table | awk '{print $2}' | sort | uniq | wc -l`
      if [[ $((N1+N2)) > $N12 ]]
      then
        echo -e "$GSE1\t$GSE2\t$((N1+N2-N12))\t$N1\t$N2"
      fi
    fi
  done
done | sort | uniq > $SPECIES.gse_xN_overlap.table

echo "Done making a list of overlapping GSEs."
awk '{if ($4==$5) {sum1++} else {sum2++}} END {printf "Same size GSEs: %d, different size GSEs: %d, overall: %d.\n",sum1,sum2,sum1+sum2}' $SPECIES.gse_xN_overlap.table

cat $SPECIES.gse_x2_overlap.table $SPECIES.gse_xN_overlap.table > $SPECIES.metrics_table2.out
