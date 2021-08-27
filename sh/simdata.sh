#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/simulations/neutpoly_sim.py

# command line args
ch=$1
init=$2
sid=$3
prm="None"

# eval optional params
for var in "${@:4}"; do 
    eval $var
done

if [ "$#" -ge 3 ] ; then
    $anaconda $prog $ch $init $sid $prm

else 
    echo "usage: neutpoly_sim <ch> <init> <sim_id> <params>"

fi

