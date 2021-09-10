#!/bin/sh

# read chrom and population IDs file from command line
CHR=$1
IDS=$2

# parse the population ID from the IDs file name
POP=$(echo $IDS | awk -F"/" '{sub(".txt", "", $NF); print $NF}')

# fix phase variable
PHASE='phase3'

# url to full VCF files
URL=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

if [ $CHR ] && [ $IDS ]; then 
    # output file path
    FOUT=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/vcf/${CHR}.${POP}.${PHASE}.BiallelicSNVs.vcf.gz
    # stream zipped VCF data from ftp URL
    curl $URL |
    # decompress VCF data and pipe to vcftools
    zcat $FIN |
    # recode biallelic snps only for the subpopulation
    vcftools --vcf - --keep $IDS --max-alleles 2 --min-alleles 2 --remove-indels --recode --recode-INFO-all --stdout |
    # compress output
    gzip --best --to-stdout > $FOUT
    
    # echo -e "chrom=$CHR\nIDs-file=$IDS\npop=$POP\nurl=$URL\nfileout=$FOUT"

else
    echo 'usage: recodesubpop <chrom> <IDs-file>'

fi

