#!/bin/bash
## Download GSE sets from NCBI using prepared file with GEO search results in summary (text) format".

if [ "$#" -ne 2 ]; then
    echo "USAGE: bash download.sh <GEO search results file path> <destination folder>"
    exit 1
fi

# path to the prepared file
METAFILE=$1

# extract URL prefixes from the file
URLS=`grep "FTP download" $METAFILE | grep GSE | perl -ne 'm/GEO( | \(.*\) )(ftp\:\/\/ftp.*\/)/; print "$2\n"' | awk -F "," '{print $1}'`

# Target directory
DIR=$2



if [ -d "$DIR" ]; then
    echo $DIR already exists
else
    mkdir -p $DIR
fi

for i in $URLS
do
    echo Processing $i
    wget -P $DIR -nd -nH -r --no-parent -A '*.gz' "$i/matrix/" 2>&1;
done
echo -e "Downloading finished"