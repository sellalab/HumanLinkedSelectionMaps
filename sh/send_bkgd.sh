#!/bin/sh
init_dir=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files

# sp=$1
# pc=$2
# pn=$3
# if [ "$#" == 3 ] ; then
#     # an=${sp}_cons${pc}_euarchontoglires_neut${pn}_filtered
#     an=${sp}_cons${pc}_Segments
#     f=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files/YRI.${sp}${pc}_euarchontoglires${pn}_filtered.BS1.4.CS0.0.NOT_STARTED.initial.txt
#     sh/bkgd_batch.sh $f n=${an} pi=0 plen=5e7 >> ${an}.bkgd.jobs.log
# else 
#     echo 'usage: send_bkgd <anno> <pcons> <pneut>'
# fi

if [ "$#" == 2 ]; then
	f=$1
	an=$2
	sh/bkgd_batch.sh $f n=${an} pi=0 plen=5e7 >> ${an}.bkgd.jobs.log
else	
	echo 'usage: send_bkgd <init> <anno>'
fi

# pc=$1
# if [ $pc ] ; then
#     # an=${sp}_cons${pc}_euarchontoglires_neut${pn}_filtered
#     an=ape_cons${pc}_euarchontoglires_neut35_filtered
#     f=${init_dir}/YRI.ape_cons${pc}_euarchontoglires_neut35_filtered.BS1.4.CS0.0.NOT_STARTED.initial.txt
#     sh/bkgd_batch.sh $f n=${an} pi=0 plen=5e7 >> ${an}.bkgd.jobs.log
# else 
#     echo 'usage: send_bkgd <pcons>'
# fi