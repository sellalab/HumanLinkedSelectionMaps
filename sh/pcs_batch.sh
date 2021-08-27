#!/bin/sh

# SEND A BATCH OF BKGD_EXACT JOBS

# path to myqsub script for submitting batch jobs
qq=~/myqsub.sh
# path to python calc_bkgd wrapper
pcs=/ifs/data/c2b2/gs_lab/dam2214/run/sh/prcs_bparts.sh
# default chrom range
ci=1
cj=22
# default tval range
# TODO: implement tval handling downstream in exact calc python
ti=2
tj=4.5
# default hours for the run
hr=1
# default memory allocation
mem=2G
# default partition length of merged file
p="1e7"
# get bdir from command line
bdir=$1

# set non-defaults from command line
for var in "${@:2}"; do 
    eval $var
done

if [ "$#" -ge 1 ] ; then 
    # print out accepted arguments
    printf "SENDING BATCH OF PCS_BATCH JOBS FOR THE FOLLOWING PARAMS:\n"
    printf "ci=%s\ncj=%s\nti=%s\ntj=%s\np=%s\nhr=%s\nmem=%s\n" $ci $cj $ti $tj $p $hr $mem

    for c in $(eval echo {$ci..$cj}) ; do 
        for t in $(seq $ti 0.5 $tj) ; do
            $qq $mem ${hr}:: pcs-batch-call.log $pcs chr${c} $bdir "t=10**-${t}" "p=${p}"
        done
    done

else
    echo -e "usage: pcs_batch <bdir> [ci=1, cj=22, ti=2, tj=4.5, p=${p} hr=${hr} mem=${mem}]";

fi






