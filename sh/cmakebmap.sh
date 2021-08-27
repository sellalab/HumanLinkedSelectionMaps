#!/bin/sh

# path to inner shell script
prog=sh/makebmap.sh

# get folder from command line
fldr=$1

mem=10G
tm=:30:

if [ $fldr ]; then 
    f_log=${fldr}.makebmap.log
    printf "%s " $fldr
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: cmakebmap <folder_name>'
fi
