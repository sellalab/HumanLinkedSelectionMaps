#!/bin/sh

# get the coverage at each neutral SNP on the chromosome
c=$1

if [ $c ]; then
    # file of neutral SNP positions
    pf=frq/YRI/neut-snps/chr${c}.neut.snp.pos.txt
    # output mean cov/SNP pos file
    cf=cov/chr${c}.neut.snp.meancov.txt

    # stream the vcxf file to stdout
    zcat vcf/chr${c}.phase-3.vcf.gz |
    # first read a SNP positions to array
    # then read complete VCF file and extract depth at neutral SNP positions in the array
    awk 'FNR==NR{a[$0]++;next}!/#/&&a[$2]{match($8,/DP=[0-9]+/);d=substr($8,RSTART+3,RLENGTH-3);n=NF-9;print $2,d/n}' $pf - > $cf

else
    echo 'usage: snpdepth <chrom>'

fi

