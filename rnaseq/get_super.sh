#!/bin/bash 

zcat *gz | grep Super | grep GSE | perl -ne 'm/(GSE.*)"/; print "$1\n"' > Super.list

KK=`cat Super.list`

for i in $KK
do
  MASK=`echo $i | perl -ne 'chomp; $l=length($_)-3; $sub=substr $_,0,$l; printf "%snnn\n",$sub'` 
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/$MASK/$i/matrix/*
done 
