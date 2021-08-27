#!/bin/sh

# path to inner shell script
prog=sh/sendbkgd.sh



mem=2G
tm=72::

if [ "$#" == 2 ]; then 
    fi=$1
    an=$2
    f_log=sendbkgd.$(basename $fi .txt).log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fi $an
elif [ "$#" == 1 ]; then 
    an=$1 
    f_log=sendbkgd.consolidate.${an}.log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $an
else
    echo 'usage_1: send_bkgd_jobs <f_init> <anno>'
    echo 'usage_2: process_partitions <anno>'
fi
