#!/usr/bin/awk -f

BEGIN{
    if (!POP){
	print "usage: allid.awk -v POP=<pop_id>"
	exit
    }
    fout=sprintf("/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps/pops/%s.all.IDs.txt", POP)
    # print fout
}
$2==POP{
    # print $1
    print $1 > fout
}
    