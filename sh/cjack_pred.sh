#!/bin/sh

# path to inner shell script
prog=sh/jack_pred.sh

# get folder from command line
anno=$1

mem=16G
tm=4::

if [ $anno ]; then 
    f_log=${anno}.jack_pred.log
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $anno
else
    echo 'usage: jackknife_prediction <anno>'
fi
