#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/compare_wigz.py

if [ "$#" == 3 ]; then 
    ch=$1 
    sp1=$2
    sp2=$3
    $py $prog $ch $sp1 $sp2
else
    echo 'usage: compare_wigz <ch> <sp1> <sp2>'
fi
