#!/bin/sh

# path to inner shell script
prog=sh/calc_rsq.sh

# get folder from command line
fldr=$1

if [ $fldr ]; then 
    f_log=${fldr}.rsq.log
    qsub -l mem=4G,time=:15: -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: crsq <folder_name>'
fi
