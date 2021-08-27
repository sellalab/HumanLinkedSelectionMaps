#!/bin/sh

# path to inner shell script
prog=sh/jack_analyses.sh

# get folder from command line
anno=$1

if [ $anno ]; then 
    f_log=${anno}.jack_analyses.log
    qsub -l mem=16G,time=4:: -cwd -j y -o $f_log $prog $anno
else
    echo 'usage: jackknife_analyses <anno>'
fi