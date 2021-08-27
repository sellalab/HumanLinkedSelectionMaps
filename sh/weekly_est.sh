#!/bin/sh

# path where gzipped accounting files are stored
acct_dir=/opt/gridengine/hpc/default/common

# input number of weeks back to be processed
wks_j=$1
# optional start week for looking at arbitrary range
wks_i=$2
# total memory counter
mem_tot=0
# free memory allocation
mem_free=$((10000*3600))

if [ $wks_i ]; then 
    if ((wks_i <= 0)); then
	echo $wks_i
	echo "ERROR: wsk_i must be >= 0"
	exit
    fi
else
    wks_i=0
fi

# print the latest usage by default unless wks_i != 0 
if [ $wks_i == 0 ] ; then 
    # get the date of the PREVIOUS save
    latest=$(ls -l $acct_dir/accounting.0.gz | awk '{print $6, $7}')

    # get memory use and estimated cost
    mem=$(qacct -o $(whoami) | awk 'NR==3{if ($5>$6) print $5; else print $6}')
    cost=$(awk -v m=$mem 'BEGIN{printf("%.2f", (5e-06 * m)); exit}')
    echo "usage since ${latest}: ${mem} GB/CPU_sec; ${cost} dollars"

    # update total memory use
    new_tot=$(awk 'BEGIN{print ARGV[1] + ARGV[2]; exit}' $mem $mem_tot)
    mem_tot=$new_tot
fi 

# print the last wks_j records if a number of weeks is given
if [ $wks_j ] ; then 
    if ((wks_j>0)) ; then 
	for i in $(eval echo {$wks_i..$wks_j}); do 
            # get the date of the PREVIOUS save
	    iplus=$((i+1))
	    latest=$(ls -l $acct_dir/accounting.${iplus}.gz | awk '{print $6, $7}')

       	    # start checking from wks_i index
	    mem=$(zcat $acct_dir/accounting.${i}.gz | grep $(whoami) | awk -F":" '{if ($37>$38) m+=$37; else m+=$38}END{print m}')
	    cost=$(awk -v m=$mem 'BEGIN{printf("%.2f", (5e-06 * m)); exit}')
	    echo "usage since ${latest}: ${mem} GB/CPU_sec; ${cost} dollars"

	    # update total memory used
	    new_tot=$(awk 'BEGIN{print ARGV[1] + ARGV[2]; exit}' $mem $mem_tot)
	    mem_tot=$new_tot

	done
	cost=$(awk -v m=$mem_tot 'BEGIN{printf("%.2f", (5e-06 * m)); exit}')
	echo "total usage for interval: ${mem_tot} GB/CPU_sec; ${cost} dollars"
    else
	echo "ERROR: wks_j must be > 0"
    fi
fi
