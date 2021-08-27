#!/bin/sh

# path to inner shell script
prog=sh/sendjk1.sh

fi=$1
istart=$2
iend=$3

mem=2G
tm=72::

if [ "$#" == 3 ]; then 
    # f_log=sendjk1.istart_${istart}_iend${iend}.log
    f_log=sendjk2.istart_${istart}_iend${iend}.log
    # f_log=cleanupjackknife.istart_${istart}_iend${iend}.log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fi $istart $iend
else
    echo 'usage: send_jackknife_jobs <f_init> <istart> <iend>'
fi
