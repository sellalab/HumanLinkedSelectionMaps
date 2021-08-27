#!/bin/sh

# path to inner shell script
prog=sh/collate_diversity.sh

if [ "$#" == 3 ]; then
    fldr=$1
    wdth=$2
    focal=$3
    f_log=${fldr}.collate.${focal}.width.${wdth}.log
    qsub -l mem=16G,time=24:: -cwd -j y -o $f_log $prog $fldr $wdth $focal

else
    echo 'usage: collate_diversity <fldr> <width> <focal>'

fi