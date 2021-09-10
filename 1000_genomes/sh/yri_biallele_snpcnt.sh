#!/bin/sh

# This script uses curl to stream zipped vcf data from ftp site. Using VCFtools the data is parsed to get SNP counts for a population subset

# remove monomorphic sites awk script
rm_mono=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/snps/sh/remove_monomorphs.awk

# read chrom and population IDs from command line
CHR=$1
POP=$2
TOKEN=$3  # all/female/male
PHASE='phase3'

if [ $CHR ] && [ $POP ] && [ $TOKEN ]; then 
    # path to individual IDs file
    # IDS=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/snps/pops/${POP}.random.81.IDs.txt
    IDS=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/snps/pops/${POP}.${TOKEN}.IDs.txt    
    # use the smaller, pre-processed YRI files for biallelic SNVs
    FIN=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/snps/vcf/yriall_biallellic/${CHR}.YRI.all.IDs.${PHASE}.BiallelicSNVs.vcf.gz
    # output file path
    FOUT=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/snps/frq/${POP}/${CHR}.${POP}.${TOKEN}.${PHASE}.frq.count    

    # pipe the unzipped VCF file to VCFtools (faster than using --gzvcf with crappy zlib)
    zcat $FIN |
    # get SNP counts with local zipped VCF file
    # vcftools --vcf - --keep $IDS --max-alleles 2 --min-alleles 2 --counts --remove-indels --stdout > $FOUT
    vcftools --vcf - --keep $IDS --max-alleles 2 --min-alleles 2 --counts --remove-indels --stdout |
    # filter monomorphic sites from stdout save to output file
    $rm_mono |
    # filter rare phase-3 duplicate SNPs using awk unique one-liner
    awk '!a[$2]++' > $FOUT

else
    echo 'snpcnt <chrom> <pop_id> <token>'

fi

