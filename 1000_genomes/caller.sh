#!/bin/sh

# path to the snpcnt script
snpcnt=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/sh/snpcnt.sh

# read chrom and population IDs from command line
CHROM=$1
POP_ID=$2

# set extra memory externally for big chroms (1-3)
if [ $3 ]; then
    MEM=$3
else
    MEM='1G'
fi

if [ $CHROM ] && [ $POP_ID ]; then 
    # run log
    log=logs/${POP_ID}/${CHROM}.${POP_ID}.snpcount.log
    # qsub command
    qsub -l mem=$MEM,time=6:: -cwd -j y -o $log $snpcnt $CHROM $POP_ID

else
    echo 'caller <chrom> <pop_id> <[mem]>'    

fi