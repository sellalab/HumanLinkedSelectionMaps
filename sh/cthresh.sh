#!/bin/sh

# path to inner shell script
prog=sh/thresh.sh

f=$1

if [ $f ]; then
	f_log=impose.thresh.log
    qsub -l mem=4G,time=:5: -cwd -j y -o $f_log $prog $f
else
    echo 'usage: impose_threshold <filename>'
fi
