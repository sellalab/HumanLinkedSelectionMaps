#!/bin/sh

# path to inner shell script
prog=sh/nxdist.sh

chrom=$1

if [ $chrom ]; then
	f_log=${chrom}.nxdist.log
    qsub -l mem=4G,time=:45: -cwd -j y -o $f_log $prog $chrom
else
    echo 'usage: nonexonic_distribution <chrom>'
fi
