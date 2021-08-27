#!/bin/sh

# path to inner shell script
prog=sh/predcor.sh

if [ "$#" == 2 ]; then
    fl1=$1
    fl2=$2
    f_log=${fl1}.${fl2}.predcor.log
    qsub -l mem=36G,time=1:: -cwd -j y -o $f_log $prog $fl1 $fl2

else
    echo 'usage: lsmap_correlation <fldr_1> <fldr_2>'

fi