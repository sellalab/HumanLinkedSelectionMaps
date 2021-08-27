#!/bin/sh

# path to inner shell script
prog=sh/jack_rsq.sh

# get folder from command line
anno=$1

mem=12G
tm=1::

if [ $anno ]; then 
    f_log=${anno}.jack_rsq.log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $anno
else
    echo 'usage: jackknife_rsq <anno>'
fi
