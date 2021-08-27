#!/bin/sh

# path to inner shell script
prog=sh/sortmore.sh

# get folder from command line
fldr=$1

mem=16G
tm=1::

if [ $fldr ]; then 
    f_log=${fldr}.sortmore.log
    printf "%s " $fldr
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: csortmore <folder_name> [mem=24G, tm=24::]'
fi
