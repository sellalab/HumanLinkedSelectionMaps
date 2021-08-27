#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/figures/nonexonic_distribution.py
chrom=$1
if [ $chrom ]; then
    $anaconda $prog $chrom
else
    echo 'usage: nonexonic_distribution <chrom>'
fi
