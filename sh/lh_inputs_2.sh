#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/lh_inputs.py

ch=$1
init=$2
plen=$3

if [[ $ch ]]; then 
    $py $prog $ch $init $plen

else
    echo 'usage: lh_inputs <chrom> <init> <plen>'

fi
