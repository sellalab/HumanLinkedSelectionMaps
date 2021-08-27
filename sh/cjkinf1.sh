#!/bin/sh

# path to inner shell script
prog=sh/jkinf1.sh

init=$1
jidx=$2

if [ "$#" == 2 ]; then
    lbl=$(grep '\-\-tkn=' $init | awk -F"=" '{gsub("\47", "", $2); print $2}')
    for idx in {0..14}; do 
    	f_log=${lbl}.idx${idx}.jkidx${jidx}.inf1.log;
		qsub -l mem=24G,time=10:: -cwd -j y -o $f_log $prog $init $idx $jidx;
    done
else
    echo "usage: cinf1 <init> <jkidx>"
fi
