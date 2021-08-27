#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/calc_rsq.py

fldr=$1

if [ $fldr ]; then 
    $py $prog $fldr
else
    echo 'usage: calc_rsq <folder_name>'

fi
