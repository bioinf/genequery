#!/bin/bash 

for i in *series_matrix.txt.gz
do
  ./parse_matrix.pl $i
done
