#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/jack_inf_step_1.py

if [ "$#" == 3 ]; then
    $anaconda $prog $1 $2 $3
else
    echo "usage: cluster_runinf <init> <idx> <jackidx>"
fi
