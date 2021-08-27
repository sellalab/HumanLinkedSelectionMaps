#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/simulations/convergence_tests.py

# get param index and method name from command line
# idx=$1
# meth=$2

init=$1
idx=$2
bth=$3

# only execute when both params are specified
if [ "$#" == 3 ]; then
    $anaconda $prog $init $idx $bth
else
    echo "usage: convergence_tests <init> <iprm> bthresh>"
fi
