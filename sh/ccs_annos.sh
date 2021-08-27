#!/bin/sh

# path to inner shell script
prog=sh/cs_annos.sh



if [ "$#" == 3 ]; then
    init=$1
    ch=$2
    expo=$3
    lbl=$(echo $init | awk -F"/" '{split($NF, a, "."); print a[2]}')	
    f_log=${ch}.${expo}.${lbl}.csmaps.log
    qsub -l mem=4G,time=16:: -cwd -j y -o $f_log $prog $init $ch $expo

elif [ "$#" == 2 ]; then
    ch=$1
    an=$2
    f_log=${ch}.${an}.HC.subs.log
    qsub -l mem=4G,time=:15: -cwd -j y -o $f_log $prog $ch $an

elif [ "$#" == 4 ]; then
    ch=$1
    an=$2
    pop=$3
    ss=$4
    f_log=${ch}.${an}.${pop}.s${ss}.HC.subs.log
    qsub -l mem=4G,time=:15: -cwd -j y -o $f_log $prog $ch $an $pop $ss

else
    echo 'usage_1: cs_annotations <chrom> <anno>'
    echo 'usage_2: cs_annotations <init> <chrom> <expo>'
    echo 'usage_3: cs_annotations <chrom> <anno> <pop> <sample>'
fi