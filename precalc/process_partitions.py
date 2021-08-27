#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 01:55:30 2018

@author: davidmurphy
"""

import os
from sys import argv
from classes.bkgdmap import BkgdMapReader
from classes.runstruct import ChromStruct
from data_processing.data_tools import f_zip


def process_partition_files(cst, anno, t, plen):
    """merge partition and save. compress partition files in subdir"""

    # create sub directory for partition files
    s_plen = '{:.0e}'.format(plen).replace('+', '')
    pdir = '{}/merge_{}'.format(cst.bdir, s_plen)
    sub_dir = '{}/precalc/{}'.format(cst.root, pdir)
    if not os.path.isdir(sub_dir):
        os.mkdir(sub_dir)

    # initialize an empty bmap reader
    bmr = BkgdMapReader()    
    for i in xrange(int(cst.chlen / plen) + 1 ):
        # original file path
        f = cst.bkgd_file(anno, t, pidx=i, plen=plen)
        
        # process with reader
        bmr.read_bfile(f)
        
        # # gzipped and moved file path
        # gz_f = f.replace(cst.bdir+'/', pdir+'/') + '.gz'
        #
        # # zip files in the sub directory (don't delete for now)
        # f_zip(f, gz_f, rm=True)

    # write final partition file
    f_out = cst.bkgd_file(anno, t, plen=plen, merge=True)
    bmr.write_bfile(f_out)


def main_remote():
    if len(argv) != 7:
        print 'usage: process_partitions <ch> <gmap> <bdir> <anno> <t> <plen>'
        exit(1)
    ch, gmap, bdir, anno = argv[1:5]
    t, plen = map(eval, argv[5:])
    cst = ChromStruct(chrom=ch, bdir=bdir, gmap=gmap)

    process_partition_files(cst, anno, t, plen)
    
        
def main_local():
    ch = 'chr2'
    bdir = 'std_split_pts'
    anno = 'primate_cons95_Segments'
    t = 10**-2
    plen = 1e6
    cst = ChromStruct(chrom=ch, bdir=bdir)
    from datetime import datetime
    for ch in 'chr1 chr2 chr3 chr8 chr9 chr10'.split():
        t_0 = datetime.now()
        cst.chrom = ch
        process_partition_files(cst, anno, t, plen)
        t_e = datetime.now() - t_0
        print '{} processed. time={}'.format(ch, t_e)


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy'):
        main_local()
    else:
        main_remote()