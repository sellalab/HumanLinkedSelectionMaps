#!/bin/sh

# path to inner shell script
prog=sh/inf1.sh

init=$1

# use default bsx and precalc bmap threshold values
min_bsx=0
min_b=0.0228202
mem=24G
tm=12::

# update variables with command line input if given
for var in "${@:2}"; do 
    eval $var
done

if [ "$#" -ge 1 ]; then
    lbl=$(echo $init | awk -F"/" '{split($NF, a, "."); print a[1]"."a[2]"."a[3]}')
    for idx in {0..14}; do 
#    for idx in 6; do 
        printf "%s.idx%s " $lbl $idx
        f_log=${lbl}.idx${idx}.min_bsx${min_bsx}.min_b${min_b}.inf1.log
        qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog $init $idx $min_bsx $min_b
    done
else
    printf 'usage: cinf1 <init> [min_bsx=%s, min_b=%s, mem=24G, tm=24::]\n' $min_bsx $min_b
fi
