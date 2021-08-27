#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/sweep_map.py

ch=$1
an=$2
sx=$3
lm=$4

if [ "$#" == 4 ]; then 
    $py $prog $ch $an $sx $lm;

else
    echo 'usage: sweep_maps chrom anno coef lim'

fi
