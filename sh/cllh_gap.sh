#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/notebook_2.py

finit=$1
bth=$2
pfld=$3


if [ $finit ]; then 
    $py $prog $finit $bth $pfld
else
    echo 'usage: notebook_2 <folder_name> [b_thresh] [param_folder]'
fi
