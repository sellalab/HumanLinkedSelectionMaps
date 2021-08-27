#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/send_jackknife_jobs.py

fi=$1
istart=$2
iend=$3

# run inf2 on folder and when that finishes, run rsq
if [ "$#" == 3 ]; then 
    $py $prog $fi $istart $iend
else
    echo 'usage: send_jackknife_jobs <f_init> <istart> <iend>'
fi
