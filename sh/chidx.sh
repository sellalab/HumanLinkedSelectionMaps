#!/bin/sh

# chromosome lengths file
lengths=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/ch_features/chr_len_all.txt

# awk script
prog=/ifs/data/c2b2/gs_lab/dam2214/run/awk/chidx.awk

# get chrom from command line
ch=$1

# get plen from command line
plen=$2

# use quiet mode to just get the number back
quiet=$3

if [ "$#" == 2 ]; then 
    # pass chrom and plen to awk script to get max idx
    awk -v ch=$1 -v plen=$2 -f $prog $lengths;

elif [ "$#" == 3 ]; then 
    awk -v ch=$1 -v plen=$2 -v quiet=1 -f $prog $lengths;
else
    echo "usage: chridx_max <ch> <plen> [q=quiet_mode]"

fi

