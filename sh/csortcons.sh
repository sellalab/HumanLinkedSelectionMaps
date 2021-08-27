#!/bin/sh

# path to inner shell script
prog=sh/sortcons.sh

# get chrom and folder from command line
ch=$1
fldr=$2
idx=$3 

mem=8G
tm=24::

if [ "$#" == 3 ]; then 
    f_log=${ch}.${fldr}.${idx}.sortcons.log
    printf "%s %s " $ch $fldr
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $ch $fldr $idx
elif [ "$#" == 1 ]; then
	f_log=${1}.sortcons.meanarr.log
	printf "%s " $1
    qsub -l mem=64G,time=6:: -cwd -j y -o $f_log $prog $1
else
 	echo 'usage_1: prediction_sorted_conservation <chrom> <folder> <idx>'
    echo 'usage_2: prediction_sorted_conservation <folder>'
fi
