#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/sample_hc_substitutions.py

if [ "$#" == 2 ]; then
    ch=$1
    pop=$2
    $anaconda $prog $ch $pop
elif [ "$#" == 4 ]; then
    ch=$1
    pop=$2
    sample=$3
    variant=$4
    $anaconda $prog $ch $pop $sample $variant 
else
    echo 'usage_1: sub_hc_substitutions <chrom> <pop>'
    echo 'usage_2: sub_hc_substitutions <chrom> <pop> <sample> <variant>'

fi