#!/bin/sh

# path to qsub script
qq=~/myqsub.sh

# path to inference run wrapper
prog=/ifs/data/c2b2/gs_lab/dam2214/run/sh/multi_runinf.sh

# path to awk script that pulls the ID from init file name
getid=/ifs/data/c2b2/gs_lab/dam2214/run/awk/get_id.awk

# optional pause step between qsub calls (seconds)
sleep_time=1

# program memory and time allocation
mem=12G
tm=3::

# init file path (required)
init=$1

# default values for other inputs
odict="None"
para="True"
new="True"
sim="False"

# ID pulled from init file used for log file name
log_id=$($getid $init)
log="runinf.${log_id}.log"

# scan command line for optional args and evaluate
echo "nondefault args:"
for var in "${@:2}"; do
	eval $var
	echo $var
done

# qsub command
if [ $init ]; then
	$qq $mem $tm $log $prog $init "odict=${odict}" "para=${para}" "new=${new}" "sim=${sim}"
	# echo -e $qq"\n"$mem"\n"$tm"\n"$log"\n"$prog"\n"$init"\n"$odict"\n"$para"\n"$new
	sleep $sleep_time


else
	echo "usage: inf_batch <init> [sleep_time=1s mem=12G time=3:: odict=None para=True new=True sim=False]"

fi

