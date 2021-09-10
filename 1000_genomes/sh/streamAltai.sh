#!/bin/sh

# awk script parses out biallelic, high quality SNVs
PROG=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/sh/biallelicHQneandertal.awk

# set chrom on command line
CHR=$1

# url for current chrom
URL=http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.${CHR}.mod.vcf.gz

# path to output file
FOUT=/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/vcf/AltaiNea.hg19_1000g.${CHR}.mod.Biallelic.HighQual.vcf

if [ $CHR ]; then
    # stream VCF file from database and parse out biallelic nonREF variants
    curl $URL | zcat | $PROG > $FOUT
else
    echo "usage: streamAltai <chrom> (numeric, e.g., '22')"
fi

