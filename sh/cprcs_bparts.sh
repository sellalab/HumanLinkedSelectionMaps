#!/bin/sh 

# program path
prog=sh/prcs_bparts.sh

ch=$1
gmap=$2
bdir=$3
anno=$4
t=$5
p=$6

if [ "$#" == 6 ]; then
	f_log=${ch}.${gmap}.${bdir}.${anno}.${t}.${p}.processbmap.log
    qsub -l mem=3G,time=:20: -cwd -j y -o $f_log $prog $ch $gmap $bdir $anno $t $p
else
    echo "usage: process_bmap <ch> <gmap> <bdir> <anno> <t> <plen>"
fi
