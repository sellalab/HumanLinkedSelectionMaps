#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/send_bkgd_jobs.py


# run inf2 on folder and when that finishes, run rsq
if [ "$#" == 2 ]; then 
    fi=$1
    an=$2
    $py $prog $fi $an
elif [ "$#" == 1 ]; then 
    an=$1
    $py $prog $an
else
    echo 'usage_1: send_bkgd_jobs <f_init> <anno>'
    echo 'usage_2: process_partitions <anno>'
fi
