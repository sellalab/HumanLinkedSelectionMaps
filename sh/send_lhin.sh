#!/bin/sh

sp=$1
pc=$2
pn=$3
pl=5e7

if [ "$#" == 3 ] ; then
    f=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files/YRI.${sp}${pc}_euarchontoglires${pn}_filtered.BS1.4.CS0.0.NOT_STARTED.initial.txt
    sh/lhin.sh $f $pl
else 
    echo 'usage: send_lhin <spec> <pcons> <pneut>'
fi

