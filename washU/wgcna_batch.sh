#!/bin/bash

## $LIST= x00..x70 - provided with qsub -v

#PBS -N cluster_$LIST
#PBS -l nodes=1:ppn=4:haswell,walltime=24:00:00
#PBS -o hs_cluster_$LIST.log
#PBS -j oe 
#PBS -q dque

module load R-3.3.2

cd /home/apredeus/GeneQuery/hs_csv

for i in `cat $LIST`
do
  TAG=${i%%_preprocessed.csv}
  echo "Processing file number $TAG.."
  while [ $(ps -A | grep -c R) -ge 4 ] ; do sleep 10; done
  Rscript /home/apredeus/GeneQuery/parsens_single_wgcna.R $i &> $TAG.cluster.log & 
done
wait
