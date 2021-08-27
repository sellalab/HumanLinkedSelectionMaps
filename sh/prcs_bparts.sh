#!/bin/sh 

# call bkgd_exact.py program

# python path
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# program path
py=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/process_partitions.py

ch=$1
gmap=$2
bdir=$3
anno=$4
t=$5
p=$6

# for var in "${@:3}"; do 
#    eval $var
# done

if [ "$#" == 6 ]; then
	# "usage: process_partitions <ch> <bdir> <anno> <t> <plen>"	
    $anaconda $py $ch $gmap $bdir $anno $t $p

else
    echo "usage: bkgd_exact <ch> <gmap> <bdir> <anno> <t> <plen>"
    
fi
