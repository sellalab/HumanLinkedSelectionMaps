#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/programs/css_exact.py

ch=$1
cd=$2
sx=$3
lm=$4

if [ "$#" == 4 ]; then 
    $py $prog $ch $cd $sx $lm;

else
    echo 'usage: css_exact ch cdir coef lim'

fi
