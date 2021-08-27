#!/bin/sh

# path to inner shell script
prog=sh/fixbmap.sh
f_bmap=$1

if [ $f_bmap ]; then
    fname=$(echo $f_bmap | awk -F"/" '{print $NF}')
    flog=${fname/.bkgd/.log}
    #echo $flog
    qsub -l mem=1G,time=:5: -cwd -j y -o $flog $prog $f_bmap
else
    echo 'usage: fix_bmaps <f_bmap>'
fi
