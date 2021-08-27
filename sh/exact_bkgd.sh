#!/bin/sh 

# call bkgd_exact.py program

# python path
anaconda=/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python

# program path
py=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/lsm_python/calibration/precision_tests.py

ch=$1
bdir=$2
anno="primate_cons95_Segments"
t="10**-4.5"
n="5000"
p="None"
rand="True"

for var in "${@:3}"; do 
    # echo $var
    eval $var
done

if [ "$#" -ge 2 ]; then
    $anaconda $py $ch $anno $bdir $t $n $p $rand

else
    # 'usage: bkgd_exact <ch> <anno> <bdir> <t> <n> <p> <rand>'
    echo "usage: bkgd_exact <ch> <bdir> [<anno=primate_cons95_Segments t=${t} n=${n} p=${p} rand=True]"
    
fi
