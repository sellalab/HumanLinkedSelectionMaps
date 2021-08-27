#!/bin/sh

# sp=$1
# pc=$2
# pn=$3

# if [ "$#" == 3 ] ; then
#     # create folder name
#     fldr=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/final_files/${sp}${pc}_euarchontoglires${pn}_filtered
#     # if it does not exist, create the directory
#     if [ ! -d $fldr ] ; then 
# 		mkdir $fldr
#     fi
#     # remove current folder contents (for reruns)
#     rm ${fldr}/*
#     # move all matching  final files to folder
#     mv /ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/final_files/YRI.${sp}${pc}_euarchontoglires${pn}_filtered* $fldr/
# else 
#     echo 'usage: make_move  <spec> <pcons> <pneut>'
# fi

tk=$1

if [ $tk ] ; then
    # create folder name
    fldr=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/final_files/${tk}
    
    # if it does not exist, create the directory
    if [ ! -d $fldr ] ; then 
		mkdir $fldr
    # remove current folder contents (for reruns)
	else
	    rm ${fldr}/*
    fi

    # move all matching final files to folder
    mv /ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/final_files/YRI.${tk}* $fldr/
    
else 
    echo 'usage: make_move <token>'
fi