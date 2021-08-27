#!/bin/sh

# call calc_bkgd via python caller
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# calc_bkgd python wrapper
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/precalc/calc_bkgd.py

# external args
i=$1
c=$2
t=$3
n=$4
x=$5
pi=$6
pl=$7
pt=$8

if [ "$#" == 8 ] ; then 
    $py $prog $i $c $t $n $x $pi $pl $pt

else
	echo "usage: calc_bkgd <init> <chrom> <coef> <anno> <fexe> <pidx> <plen> <pts>"
    # echo "usage: calc_bkgd <init> <chrom> <tval> <anno> <fexe> <pts>"

fi
