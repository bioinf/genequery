#!/bin/bash
## Download GSE sets from NCBI using prepared file with GEO search results in summary (text) format".
## ./download_gse.sh <GEO search results file path>"

# path to the prepared file
FILE=$1

# extract URL prefixes from the file
URLS=`grep "FTP download" $FILE | grep GSE | perl -ne 'm/GEO( | \(.*\) )(ftp\:\/\/ftp.*\/)/; print "$2\n"' | awk -F "," '{print $1}'`

for i in $URLS
do
    echo "Downloading from address $i" 
    wget $i/matrix/*_series_matrix.txt.gz  
done


echo "All downloads are now complete!" 
