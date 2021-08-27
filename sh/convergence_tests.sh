#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/simulations/convergence_tests.py

# get param index and method name from command line
# idx=$1
# meth=$2

init=$1
m=$2
d=$3
s=$4
i=$5
r=$6
f=$7
a=$8
bth=$9

# only execute when both params are specified
if [ "$#" == 9 ]; then
    $anaconda $prog $init $m $d $s $i $r $f $a $t $bth
else
    echo "usage: convergence_tests <init> <meth> <dist> <sample> <iprm> <rnd> <fudel> <alpha> <bthresh>"
fi
