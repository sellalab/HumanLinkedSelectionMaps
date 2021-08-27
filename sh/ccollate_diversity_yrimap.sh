#!/bin/sh

# path to inner shell script
prog=sh/collate_diversity_yrimap.sh

if [ "$#" == 4 ]; then
    fldr=$1
    wdth=$2
    focal=$3
    mpop=$4
    f_log=${fldr}.collate.${focal}.width.${wdth}.mpop${mpop}.log
    qsub -l mem=16G,time=24:: -cwd -j y -o $f_log $prog $fldr $wdth $focal $mpop

else
    echo 'usage: collate_diversity <fldr> <width> <focal> <map_pop>'

fi