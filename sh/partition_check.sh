#!/bin/sh

s=std_split/AA_Map_10kb.chr${c}.primate_cons95_Segments.t0.01000000.p*.1e07.bkgd
for c in {1..22}; do 
    for f in $s; do 
	awk '!/#/{n+=$2}END{print n}' $f ; 
    done | awk -v ch=chr${c} '{a[$0]++}END{print ch; for (i in a) print i, a[i]; print "----"}' ;
done