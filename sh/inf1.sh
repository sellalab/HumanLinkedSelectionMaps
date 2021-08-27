#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/inf_step_1.py

# # use default bsx and precalc bmap threshold values
# min_bsx=0.65
# min_b=None

# # update variables with command line input if given
# for var in "${@:3}"; do 
# 	eval $var
# done

if [ "$#" == 4 ]; then
    $anaconda $prog $1 $2 $3 $4
else
    echo 'usage: cluster_runinf <init> <idx> <min_bsx> <min_b'
fi
