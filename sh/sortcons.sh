#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/prediction_sorted_conservation.py

ch=$1
fldr=$2
idx=$3 

# run inf2 on folder and when that finishes, run rsq
if [ "$#" == 3 ]; then 
    $py $prog $ch $fldr $idx
elif [ "$#" == 1 ]; then
	$py $prog $1
else
    echo 'usage_1: prediction_sorted_conservation <chrom> <folder> <idx>'
    echo 'usage_2: prediction_sorted_conservation <folder>'
fi
