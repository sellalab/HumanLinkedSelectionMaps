#!/bin/sh

# path to inner shell script
prog=sh/yrisortpred.sh

# get folder from command line
fldr=$1

mem=16G
tm=:30:

if [ $fldr ]; then 
    f_log=${fldr}.yrisortpred.log
#    printf "%s " $fldr
    qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $fldr
else
    echo 'usage: sort_pred_yri <folder_name>' 
fi
