#!/bin/sh

# SEND A BATCH OF BKGD_EXACT JOBS

# path to myqsub script for submitting batch jobs
qq=~/myqsub.sh
# path to python calc_bkgd wrapper
xct=/ifs/data/c2b2/gs_lab/dam2214/run/sh/exact_bkgd.sh
# default chrom range
ci=1
cj=22
# default tval range
# TODO: implement tval handling downstream in exact calc python
ti=4.5
tj=4.5
# default hours for the run
hr=1
# default memory allocation
mem=2G
# default sample size for random points
n=5000
# default partition length of merged file
p=None
# get bdir from command line
bdir=$1

# set non-defaults from command line
for var in "${@:2}"; do 
    eval $var
done

if [ "$#" -ge 1 ] ; then 
    # print out accepted arguments
    printf "SENDING BATCH OF BKGD_EXACT JOBS FOR THE FOLLOWING PARAMS:\n"
    printf "ci=%s\ncj=%s\nti=%s\ntj=%s\nn=%s\np=%s\nhr=%s\nmem=%s\n" $ci $cj $ti $tj $n $p $hr $mem

    for c in $(eval echo {$ci..$cj}) ; do 
        for t in $(seq $ti 0.5 $tj) ; do
            $qq $mem ${hr}:: xct-batch-call.log $xct chr${c} $bdir "t=10**-${t}" "n=${n}" "p=${p}"
            # echo $xct chr${c} $bdir "tval=10**-${t}" "n=${n}"
        done
    done

else
    echo -e "usage: xct_batch <bdir> [ci=1, cj=22, ti=4.5, tj=4.5, n=${n} p=${p} hr=${hr} mem=${mem}]";

fi






