#!/bin/sh

# zip a snpcount file

f=$1

if [ $f ]; then 
    gzip -9 $f

else
    echo 'usage: zipsnpcount <snpcount_file>'

fi
