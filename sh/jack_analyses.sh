#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# paths to python files
prog1=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/jackknife_prediction.py
prog2=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/jackknife_params.py
prog3=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/jackknife_rsq.py

anno=$1

if [ $anno ]; then 
    # get params
    $py $prog1 $anno
    # get leave one out prediction map
    $py $prog2 $anno
    # get rsq from prediction map
    $py $prog3 $anno

else
    echo 'usage: jackknife_analyses <anno>'

fi
