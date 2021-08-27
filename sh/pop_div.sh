#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/pop_div_sort.py

if [ "$#" == 2 ]; then 
    $py $prog $1 $2
else
    echo 'usage: pop_div_sort <pop_1> <pop_2>'
fi
