#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/impose_threshold.py
f=$1
if [ $f ]; then
    $anaconda $prog $f
else
    echo 'usage: impose_threshold <filename>'
fi
