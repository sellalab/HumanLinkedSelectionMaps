#!/bin/sh

# SCRIPT PATHS
# path to myqsub script for submitting batch jobs
qq=~/myqsub.sh
# path to python calc_bkgd wrapper
bs=/ifs/data/c2b2/gs_lab/dam2214/run/sh/calc_bkgd.sh
# path to executable in C
fx=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/debug_bkgd/src/calc_bkgd

# FIXED PARAMS
# default enforced points setting
pts=0
# default partition length
plen=5e7
# default hours for the run
hr=24
# default memory allocation
mem=1G

# VARIABLE PARAMS
ch=$1
sp=$2
pc=$3
pn=$4
ti=$5
pi=$6
an=${sp}_cons${pc}_euarchontoglires_neut${pn}_filtered
f=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files/YRI.${sp}${pc}_euarchontoglires${pn}_filtered.BS1.4.CS0.0.NOT_STARTED.initial.txt
	
# init file
# f=$1
# chrom
# ch=$2
# annotation
# an=$3
# coeff exponent
# ti=$4
# block index
# pi=$5

# rerun the failed runs
if [ "$#" == 6 ]; then
	$qq $mem ${hr}:: bkgd-rerun-call.log $bs $f $ch 10**-${ti} $an $fx $pi $plen $pts
else
	echo "bkgd_rerun <chrom> <spec> <pcons> <pneut> <tval> <pidx>"
fi
