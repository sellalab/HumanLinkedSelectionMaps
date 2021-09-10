#!/usr/bin/awk -f

# print header
/^#/{
    print
}
# print biallelic nonREF SNV that are NOT LowQual flagged
$0!~/^#/ && $4~/^[ACGT]$/ && $5~/^[ACGT]$/ && $7 !~ /LowQual/{
    print
}