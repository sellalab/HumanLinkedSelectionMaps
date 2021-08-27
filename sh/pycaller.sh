#!/bin/sh

# a multi-purpose shell script to wrap python code for qsub
# replaces individual scripts used for each python program

# python path
ANACONDA=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# the first position should be the path to a python program
PROGRAM=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/caller.py
# PROGRAM=$1

# the program is then executed including any additional arguments
$ANACONDA $PROGRAM "${@:1}"

