#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to python file
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/jackknife_params.py

anno=$1

if [ $anno ]; then 
    $py $prog $anno		
else
    echo 'usage: jackknife_params <anno>'
fi
