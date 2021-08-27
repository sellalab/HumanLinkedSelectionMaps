#!/bin/sh

sp=$1
pc=$2
pn=$3
bt=0.65

if [ "$#" == 3 ] ; then
    f=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files/YRI.${sp}${pc}_euarchontoglires${pn}_filtered.BS1.4.CS0.0.NOT_STARTED.initial.txt
    sh/cinf1.sh $f $bt
else 
    echo 'usage: send_cinf1 <spec> <pcons> <pneut>'
fi

