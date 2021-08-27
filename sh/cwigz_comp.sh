#!/bin/sh

# program path
prog=sh/wigz_comp.sh

if [ "$#" == 3 ]; then 
    ch=$1 
    sp1=$2
    sp2=$3
    # f_log=${ch}.${sp1}.${sp2}.wigz_comp.log
    f_log=${ch}.${sp1}.${sp2}.measure_missing.log

    qsub -l mem=8G,time=2:: -cwd -j y -o $f_log $prog $ch $sp1 $sp2
else
    echo 'usage: compare_wigz <ch> <sp1> <sp2>'
fi
