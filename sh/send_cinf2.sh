#!/bin/sh

sp=$1
pc=$2
pn=$3

if [ "$#" == 3 ] ; then
    fldr=${sp}${pc}_euarchontoglires${pn}_filtered
    sh/cinf2.sh $fldr
else 
    echo 'usage: send_cinf2 <spec> <pcons> <pneut>'
fi

