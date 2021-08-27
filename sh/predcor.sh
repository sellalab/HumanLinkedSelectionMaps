#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/lsmap_correlation.py

if [ "$#" == 2 ]; then
	fl1=$1
	fl2=$2
    $anaconda $prog $fl1 $fl2

else
	echo 'usage: lsmap_correlation <fldr_1> <fldr_2>'

fi