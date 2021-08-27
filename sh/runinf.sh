#!/bin/sh 

# call run_inference.py program

# python path
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# program path
py=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/run_inference.py

bdir=$1
scale="100.0"

# can set scale with kwd
for var in "${@:2}"; do 
    eval $var
done

if [ "$#" -ge 1 ]; then
    $anaconda $py $bdir $scale

else
    echo "usage: run_inference <bdir> <bscale>"
    
fi
