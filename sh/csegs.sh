#!/bin/sh

# path to inner shell script
prog=sh/segs.sh

fldr=$1

if [ $fldr ]; then
    # f_log=${ch}.${init}.segdata.log
    for c in {1..22}; do     
    	f_log=chr${c}.${fldr}.segdata.log
	qsub -l mem=12G,time=2:: -cwd -j y -o $f_log $prog $fldr chr${c}
    done	
else
    echo 'usage: segmented_data <folder>'
fi
