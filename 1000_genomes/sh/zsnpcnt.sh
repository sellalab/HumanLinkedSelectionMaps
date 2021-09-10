#!/bin/sh

# This script uses curl to stream zipped vcf data from ftp site. Using VCFtools the data is parsed to get SNP counts for a population subset

# read chrom and population IDs from command line
CHROM=$1
POP_ID=$2

if [ $CHROM ] && [ $POP_ID ]; then 
    # path to individual IDs file
    SUB_POP=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/pops/${POP_ID}.random.81.IDs.txt
    # output file path
    OUT_FILE=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/frq/${POP_ID}/${CHROM}.${POP_ID}.frq.count
    # input file path /ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/vcf/chr11.phase-3.vcf.gz
    IN_FILE=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/vcf/${CHROM}.phase-3.vcf.gz
    
    # get SNP counts with from zipped VCF
    vcftools --gzvcf $IN_FILE --keep $SUB_POP --max-alleles 2 --min-alleles 2 --counts --remove-indels --stdout |
    # filter monomorphic sites from stdout (where either REF or ALT is equal to the sample size sample size) - save to output file
    ./remove_monomorphs.awk > $OUT_FILE

else
    echo 'snpcnt <chrom> <pop_id>'

fi

