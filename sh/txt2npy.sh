#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/txt2npy.py

if [ $1 ]; then 
    $py $prog $1

else
    echo 'usage: txt2npy <txt_file>'

fi
