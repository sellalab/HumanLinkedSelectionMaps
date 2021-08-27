#!/bin/sh

# full path to anaconda: 
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sort_predictions.py
prog=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/data_processing/sample_conserved.py

# if [ "$#" == 2 ]; then
#     $anaconda $prog $1 $2
# else
#     # echo 'usage: sample_conserved <chrom> <cons> <pct> <npct>'
#     echo 'usage: sample_conserved <ch> <cons>'    
#     # echo 'usage: sample_conserved <folder_name>'
# fi

if [ "$#" == 4 ]; then
    $anaconda $prog $1 $2 $3 $4
else
    echo 'usage: sample_conserved <chrom> <cons> <pmin> <pmax>'
fi

# if [ "$#" == 2 ]; then
#     $anaconda $prog $1 $2
# else
#     echo 'usage: sample_conserved <chrom> <cons>'
# fi

# if [ "$#" == 3 ]; then
#     $anaconda $prog $1 $2 $3
# else
#     echo 'usage: sample_conserved <chrom> <pct> <bin_size>'
# fi