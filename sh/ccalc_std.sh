#!/bin/sh

# path to inner shell script
prog=sh/calc_std.sh

# get folder from command line
fldr=$1

mem=12G
tm=:45:


if [ $fldr ]; then 
    f_log=${fldr}.std.log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: calc_std <folder_name>'
fi
