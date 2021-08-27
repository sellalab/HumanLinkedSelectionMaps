#!/bin/sh

# change dirs to the work dir
cd /ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_matlab/run

# path to matlab binary files
mat=/nfs/apps/matlab/2016a/bin/matlab

# default input script to call with matlab
prg=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_matlab/run/runInf_template.m

# matlab run options 
opt="-nojvm -nodisplay -nodesktop -singleCompThread"

# interpret any command line args...
# for var in "${@:1}"; do 
#     eval $var
# done

# run the matlab script
$mat $opt < $prg
# $mat $opt < ${prg}
# $mat -nojvm -nodisplay -nosplash -singleCompThread < ${prg}

