#!/bin/sh

# path to inner shell script
prog=sh/samplesubs.sh

if [ "$#" == 2 ]; then
    ch=$1
    pop=$2
    f_log=${ch}.${pop}.samplesubs.log
    qsub -l mem=4G,time=:45: -cwd -j y -o $f_log $prog $ch $pop

elif [ "$#" == 4 ]; then
    ch=$1
    pop=$2
    sample=$3
    variant=$4
    f_log=${ch}.${pop}.${variant}subs.sample.${sample}.log
    qsub -l mem=16G,time=:45: -cwd -j y -o $f_log $prog $ch $pop $sample $variant

else
    echo 'usage_1: sub_hc_substitutions <chrom> <pop>'
    echo 'usage_2: sub_hc_substitutions <chrom> <pop> <sample> <variant>'

fi