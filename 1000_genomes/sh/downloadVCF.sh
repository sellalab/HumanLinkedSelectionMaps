#!/bin/sh

ch=$1
if [ $ch ]; then
	
	if [[ $ch = "chrX" ]]; then 
		url=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
	
	else
		v=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		url=${v/chr1/$ch}
	
	fi	

    curl  $url > ${ch}.phase3.vcf.gz

else
    echo 'usage: downloadVCF <ch>'

fi