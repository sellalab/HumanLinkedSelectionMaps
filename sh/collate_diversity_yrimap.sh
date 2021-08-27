#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/results/collate_diversity.py

fldr=$1
wdth=$2
focal=$3
mpop=$4

if [ "$#" == 4 ]; then 
    $py $prog $fldr $wdth $focal $mpop
else
    echo 'usage: collate_diversity <folder_name> <width> <focal> <map_pop>'

fi
