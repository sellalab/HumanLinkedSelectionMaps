#!/bin/sh

# program path
prog=sh/pop_div.sh

if [ "$#" == 2 ]; then 
    f_log=${1}.${2}.popdivsort.log
    qsub -l mem=8G,time=2:: -cwd -j y -o $f_log $prog $1 $2
else
    echo 'usage: pop_div_sort <pop_1> <pop_2>'
fi
