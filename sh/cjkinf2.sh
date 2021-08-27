#!/bin/sh

# path to inner shell script
prog=sh/jkinf2.sh

# get folder from command line
fldr=$1

mem=24G
tm=10::

for var in "${@:2}"; do 
    eval $var
done

if [ $fldr ]; then 
    f_log=${fldr}.inf2.log
    printf "%s " $fldr
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: cjkinf2 <folder_name> [mem=24G, tm=24::]'
fi
