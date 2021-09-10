#!/bin/sh

# remove duplicate SNPs if they occurr

snpfile=$1
tmpfile=${snpfile}.tmp

if [ $snpfile ]; then
    # awk unique check to temp file
    awk '!a[$2]++' $snpfile > $tmpfile
    # overwrite snp file with unique sites file
    mv $tmpfile $snpfile
    
else
    echo 'usage: rmdups <SNPfile>'
fi

