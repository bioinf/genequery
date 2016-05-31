#!/bin/bash

# Download GPL annotation from NCBI.
# Example of use (on 4 processors): < GPL.list | awk '{ print $2 }' | xargs -P 4 -I{} sh -c 'bash download_gpl_annot.sh {} GPLs'
# where GPLs - folder download in, GPL.list - file of form "<n> <GPL>".

GPL=$1
DIR=$2

# GPL1234 => ../GPL1nnn/..
# GPL5432 => ../GPL5nnn/..
# GPL###  => ../GPLnnn/..
if [ ${#GPL} -eq 8 ]; then
    URL="ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL${GPL:3:2}nnn/$GPL/annot/$GPL.annot.gz"
elif [ ${#GPL} -eq 7 ]; then
    URL="ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL${GPL:3:1}nnn/$GPL/annot/$GPL.annot.gz"
else
    URL="ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/$GPL/annot/$GPL.annot.gz"
fi

wget -O $DIR/$GPL.annot.gz $URL

