#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/cs_annotations.py

if [ "$#" == 3 ]; then
	init=$1
	ch=$2
	expo=$3
    $anaconda $prog $init $ch $expo

elif [ "$#" == 2 ]; then
	ch=$1
	an=$2
    $anaconda $prog $ch $an
elif [ "$#" == 4 ]; then
	ch=$1
	an=$2
	pop=$3
	ss=$4
	$anaconda $prog $ch $an $pop $ss
else
	echo 'usage_1: cs_annotations <chrom> <anno>'
    echo 'usage_2: cs_annotations <init> <chrom> <expo>'
    echo 'usage_3: cs_annotations <chrom> <anno> <pop> <sample>'
fi