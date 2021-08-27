#!/bin/sh

# path to inner shell script
prog=sh/samplecons.sh

# ch=$1
# cons=$2
# pct=$3
# npct=$4
# fldr=$1

# if [ "$#" == 2 ]; then
# 	# f_log=${fldr}.checklengths.log
#  #    qsub -l mem=2G,time=:30: -cwd -j y -o $f_log $prog $fldr
#     f_log=${ch}.${cons}.filtered.dist.log
#     qsub -l mem=8G,time=:30: -cwd -j y -o $f_log $prog $ch $cons
#     # f_log=${ch}.${cons}${pct}.euar${npct}.create_segs.log
#     # qsub -l mem=8G,time=:30: -cwd -j y -o $f_log $prog $ch $cons $pct $npct
# else
#     echo 'usage: sample_conserved <ch> <cons>'
# fi

ch=$1
cons=$2
pmin=$3
pmax=$4
if [ "$#" == 4 ]; then
    f_log=${ch}.${cons}${pmin}to${pmax}.euar35.create_segs.log
    qsub -l mem=8G,time=:30: -cwd -j y -o $f_log $prog $ch $cons $pmin $pmax
else
    echo 'usage: sample_conserved <chrom> <cons> <pmin> <pmax>'
fi

# ch=$1
# cons=$2
# if [ "$#" == 2 ]; then
#     f_log=${ch}.${cons}.fish.zeros.dist.log
#     qsub -l mem=8G,time=:30: -cwd -j y -o $f_log $prog $ch $cons
# else
#     echo 'usage: sample_conserved <chrom> <cons>'
# fi

# ch=$1
# pct=$2
# bin_size=$3
# if [ "$#" == 3 ]; then
#     f_log=${ch}.${pct}.gmap.dist.bin_${bin_size}.log
#     qsub -l mem=12G,time=:30: -cwd -j y -o $f_log $prog $ch $pct $bin_size
# else
#     echo 'usage: sample_conserved <chrom> <pct> <bin_size>'
# fi