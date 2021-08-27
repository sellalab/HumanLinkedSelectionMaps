#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/xaut_ratio.py

# command line args
ch=$1
nt=$2 

if [ $ch ] && [ $nt ]; then
    $anaconda $prog $ch $nt

else
    echo "usage: xaut_masks <ch> <neut>"

fi
