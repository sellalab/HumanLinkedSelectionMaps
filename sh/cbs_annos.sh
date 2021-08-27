#!/bin/sh

# path to inner shell script
prog=sh/bs_annos.sh

ch=$1
pct=$2
if [ "$#" == 2 ]; then
    f_log=${ch}.${pct}.exnex.cons.log
    qsub -l mem=4G,time=:30: -cwd -j y -o $f_log $prog $ch $pct
else
    echo 'usage: bs_annotations <chrom> <anno>'
fi