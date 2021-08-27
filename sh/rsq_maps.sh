#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to compress_neutrals.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/results/rsq_maps.py

# command line args
init=$1

# new="True"
# sim="False"

# for var in "${@:2}" ; do 
#    eval $var
# done

if [ $init ]; then 
#    $anaconda $prog $init $new $sim
  $anaconda $prog $init
else
#    echo 'usage: rsquared <init> <new> <sim>'
    echo 'usage: rsquared <init>'
fi
