#!/bin/sh

# path to inner shell script
prog=sh/lh_inputs_2.sh

# fix plen for now
plen=2e7

# get the other args from command line
# ch=$1
init=$1
plen=$2
mem=16G

if [ "$#" == 2 ]; then
    for c in {1..22}; do
	ch=chr${c}
    lbl=$(echo $init | awk -F"/" '{split($NF, a, "."); print a[1]"."a[2]}')
	# lbl=$(grep '\-\-tkn=' $init | awk -F"=" '{gsub("\47", "", $2); print $2}')
	f_log=${ch}.${lbl}.lhin.log
	qsub -l mem=${mem},time=1:: -cwd -j y -o $f_log $prog $ch $init $plen
    done
else
    echo 'usage: lh_inputs <init> <plen>'
fi