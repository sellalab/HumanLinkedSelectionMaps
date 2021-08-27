#!/bin/sh

# path to inner shell script
prog=sh/nebt.sh

fldr=$1

if [ "$#" -ge 1 ]; then
    f_log=${fldr}.nebt.log
    qsub -l mem=4G,time=:15: -cwd -j y -o $f_log $prog $fldr
    
else
    echo 'usage: nebt <folder>'
fi
