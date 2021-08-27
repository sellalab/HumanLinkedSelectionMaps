#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/runinf.py

# command line args
init=$1

# optimization "options" dict string (optional)
odict="None"
para="True"
new="True"
sim="False"

# eval optional params
for var in "${@:2}"; do 
    eval $var
done


if [ "$#" -ge 1 ] ; then 
    # echo $odict
    $anaconda $prog $init $odict $para $new $sim

else
    echo "usage: runinf <init> <odict=None> <parallel=True> <new=True> <sim=False>"

fi
