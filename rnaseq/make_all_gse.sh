#!/bin/bash


KK=`cat GSE.list`
COUNT=0

for i in $KK
do
  while [ $(jobs | wc -l) -ge 16 ] ; do sleep 5; done
  ./make_gse_matrix.pl $i Rat_triple_table.tsv rnor_v6.gene_tr.tsv rnor_v6.gene_name_entrez.tsv Rat_GSE_selected_v2.tsv > $i.counts.tsv & 
  #echo "Processing matrix $COUNT, series ID $i"
  COUNT=$((COUNT+1))
done 

wait  

