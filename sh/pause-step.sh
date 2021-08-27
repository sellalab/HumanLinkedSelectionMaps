#!/bin/sh

jcount=$(qstat|wc -l)
while [ $jcount -ge 900 ]; do 
	echo "job count = ${jcount}. pausing 5 min."
	sleep 2m
	jcount=$(qstat|wc -l)
done	

