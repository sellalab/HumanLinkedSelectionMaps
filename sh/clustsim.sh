#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/simulations/cluster_simulations.py

# command line args
init=$1
ch=$2
dst=$3
smp=$4
rnd=$5
fu=$6
alp=$7
bthr=$8

if [ "$#" == 8 ] ; then
    $anaconda $prog $init $ch $dst $smp $rnd $fu $alp $bthr

else
    echo "usage: cluster_simulations <init> <ch> <dist> <sample> <random> <f_udel> <alpha> <bthresh>"

fi
