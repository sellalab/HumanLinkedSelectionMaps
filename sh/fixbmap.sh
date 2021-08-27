#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to python script
# prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/fix_bmaps.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/fill_tail_gaps.py


if [ $1 ]; then
    $anaconda $prog $1
else
    echo 'usage: fix_bmaps <f_bmap>'
fi
