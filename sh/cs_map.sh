#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/classes/cscalculator.py

ch=$1
sx=$2

if [ "$#" == 2 ]; then 
    $py $prog $ch $sx;

else
    echo 'usage: cscalculator <chrom> <coef>'

fi
