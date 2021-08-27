#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/neutral_data.py

# run python program
if [ "$#" == 3 ]; then
    # chrom outg type (b, cm)
    $anaconda $prog $1 $2 $3
elif [ "$#" == 4 ]; then
    # chrom neut outg tkn
    $anaconda $prog $1 $2 $3 $4
elif [ "$#" == 8 ]; then
    # neut outg tkn pi div b dist
    $anaconda $prog $1 $2 $3 $4 $5 $6 $7 $8
else
    echo "usage_1: get_arrays <chrom> <outg> <type>"
    echo "usage_2: get_arrays <chrom> <neut> <outg> <tkn>"
    echo "usage_3: join_arrays <neut> <outg> <tkn> <x> <pi> <div> <b> <dist>"
fi

