#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog_2=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/calc_rsq_fixmap.py

fldr=$1

# run inf2 on folder and when that finishes, run rsq
if [ $fldr ]; then 
    $py $prog_2 $fldr
else
    echo 'usage: mean_best <folder_name>'
fi
