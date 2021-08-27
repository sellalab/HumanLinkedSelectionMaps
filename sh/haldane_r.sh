#!/bin/sh 

# call bkgd_exact.py program

# python path
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# program path
py=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/programs/haldane_r.py

ch=$1
cf=$2
bd=$3


if [ "$#" -eq 3 ]; then
    $anaconda $py $ch $cf $bd

else
    echo 'usage: haldane_r ch coeff bdir'
    
fi

