#!/bin/sh

# path to inner shell script
prog=sh/jack_prms.sh

# get folder from command line
anno=$1

mem=2G
tm=1::

if [ $anno ]; then 
    f_log=${anno}.jack_prms.log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $anno
else
    echo 'usage: jackknife_params <anno>'
fi
