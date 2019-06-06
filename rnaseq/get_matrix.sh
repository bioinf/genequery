#!/bin/bash 

TXT=$1

grep FTP $TXT | perl -ne 'm/(ftp.*)/; print "$1\n"' > ftp.list

KK=`cat ftp.list`

for i in $KK
do
  ## so apparently this shit works 
  wget $i/matrix/*
done 
