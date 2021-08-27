#!/bin/sh

# sp=$1
# pc=$2
# pn=$3
# wn=5e7
# # bkgd_exact <ch> <bdir> <anno> <t> <plen>
# if [ "$#" == 3 ] ; then
#     for t in 2 2.5 3 3.5; do
# 	for c in {1..22}; do
# 	    ch=chr${c}
# 	    an=${sp}_cons${pc}_euarchontoglires_neut${pn}_filtered
# 	    qsub -l mem=1G,time=:5: -cwd -j y -o ${ch}.${an}.${t}.processbmap.log sh/prcs_bparts.sh $ch $an $an 10**-${t} $wn > /dev/null
# 	done;
#     done;
# else
#     echo 'usage: process_bmaps <spec> <pcons> <pneut>'
# fi

an=$1
bd=$2
wn=5e7
# bkgd_exact <ch> <bdir> <anno> <t> <plen>
if [ "$#" == 2 ] ; then
    for t in 2 2.5 3 3.5; do
		for c in {1..22}; do
		    ch=chr${c}
		    qsub -l mem=1G,time=:5: -cwd -j y -o ${ch}.${an}.${t}.processbmap.log sh/prcs_bparts.sh $ch $bd $an 10**-${t} $wn > /dev/null
		done;
    done;
else
    echo 'usage: process_bmaps <anno> <bdir>'
fi