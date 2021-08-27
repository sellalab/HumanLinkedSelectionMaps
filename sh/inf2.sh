#!/bin/sh

# full path to anaconda: 
py=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# full path to sweep_map.py
prog_1=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/inf_step_2.py
prog_2=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/calc_rsq.py
prog_3=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/sort_pred.py
prog_4=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/chrom_map.py
prog_5=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/results/collate_diversity.py
prog_6=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/likelihood/construct_bmap.py

fldr=$1

# run inf2 on folder and when that finishes, run rsq
if [ $fldr ]; then 
    $py $prog_1 $fldr
    $py $prog_2 $fldr
    $py $prog_3 $fldr
    $py $prog_4 $fldr
    $py $prog_5 $fldr 5e-3 nonsyn
    $py $prog_6 $fldr
else
    echo 'usage: mean_best <folder_name>'

fi
