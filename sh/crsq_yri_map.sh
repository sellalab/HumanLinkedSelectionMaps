#!/bin/sh

# path to inner shell script
prog=sh/rsq_yri_map.sh

# get folder from command line
fldr=$1

mem=24G
tm=24::

for var in "${@:2}"; do 
    eval $var
done

if [ $fldr ]; then 
    f_log=${fldr}.inf2.log
    printf "%s " $fldr
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: crsq_yri_map <folder_name> [mem=24G, tm=24::]'
fi
