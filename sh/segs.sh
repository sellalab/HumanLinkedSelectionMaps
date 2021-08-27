#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/segmented_data.py

if [ "$#" == 2 ]; then
    $anaconda $prog $1 $2    
else
    echo 'usage: segmented_data <folder> <chrom>'
fi

