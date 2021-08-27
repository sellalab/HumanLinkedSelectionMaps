#!/bin/sh

# SEND A BATCH OF CALC_BKGD JOBS BY SETTING TVALS AND CHROM-RANGE

# path to myqsub script for submitting batch jobs
qq=~/myqsub.sh

# path to python calc_bkgd wrapper
bs=/ifs/data/c2b2/gs_lab/dam2214/run/sh/calc_bkgd.sh

# path to index file: retrieves chrom idx count given plen
chidx=/ifs/data/c2b2/gs_lab/dam2214/run/sh/chidx.sh

# required external vars (init file only)
f=$1

# default anno
n='primate_cons95_Segments'

# default chrom range
ci=1
cj=22

# default tval range
ti=2
tj=4.5

# default pidx range (triggers default behavior without indexing)
pi=-1
pj=-1

# default executable file to use to create bmap
# fx=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/bkgd/src/calc_bkgd
# use the debugging bkgd as default (print warning anytime this is in use)
fx=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/debug_bkgd/src/calc_bkgd
echo "WARNING! THE DEFAULT CALC_BKGD EXECUTABLE IS THE DEBUGGING VERSION!"

# default enforced points setting
pts=0

# default partition length
plen=1e7

# default hours for the run
hr=24

# default memory allocation
mem=1G

# default interval for tvals
ivl=0.5

# set non-defaults from command line
for var in "${@:2}"; do 
    if [[ $var == fx* ]]; then
	echo "---> CALC_BKGD EXECUTABLE WAS UPDATED"
    fi
    eval $var
done

if [ "$#" -ge 1 ] ; then 
    # print out accepted arguments
    printf "SENDING BATCH OF CALC_BKGD JOBS FOR THE FOLLOWING PARAMS:\n"
    sfmt="init=%s\nanno=%s\nci=%s\ncj=%s\nti=%s\ntj=%s\npi=%s\npj=%s\nplen=%s\nhr=%s\nfx=%s\npts=%s\nmem=%s\nivl=%s\n"
    printf $sfmt $f $n $ci $cj $ti $tj $pi $pj $plen $hr $fx $pts $mem $ivl

    for c in $(eval echo {$ci..$cj}) ; do 
        # get the chromosome-specific max pidx by default if pi is specified and pj is not
        if [ $pi -ge 0 ] && [ $pj == -1 ]; then
            cmd="${chidx} chr${c} ${plen} -q"
            pjdx=$($cmd)
        else
            pjdx=$pj
        fi
    	for t in $(seq $ti $ivl $tj) ; do 
            for p in $(seq $pi 1 $pjdx) ; do 
            	# print string to appear behind job number with relevant params for tracing failed jobs
            	printf "%s %s %s: " chr${c} $t $p 
        	# wrap call with qsub caller script to submit parallel jobs
        	$qq $mem ${hr}:: bkgd-batch-call.log $bs $f chr${c} 10**-${t} $n $fx $p $plen $pts #; sleep 0.5  #| awk '{for (i=1;i<=NF;i++) print $i}'
            done
    	done
    done

else
    echo -e "usage: bkgd_batch <init> [n=anno, ti=2.0, tj=4.5, ci=1, cj=22, pi=-1, pj=-1, plen=1e7, hr=24, fx=fexe, pts=0/1, mem=mem, ivl=ivl]";

fi

