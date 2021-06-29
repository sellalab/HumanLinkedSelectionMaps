__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from itertools import izip
from data_processing.functions import swap_root
from classes.runstruct import ChromStruct, root_dir


def main():
    if len(argv) != 2:
        print 'usage: pop_div_sort <folder_name>'
        exit(1)

    fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # create a bmap dir in the result dir if it does not exist
    bdir = fdir + 'bmap/'
    if not os.path.isdir(bdir):
        os.mkdir(bdir)

    # create bmap for each chromosome
    for ch in cst.chroms:
        # reset chrom
        cst.chrom = ch
        # load segment lengths for chrom
        sg = np.load(cst.sg_files)['sg'].astype(int)
        assert np.sum(sg) == cst.chlen
        # load bmap matrix for chrom
        bs = np.load(cst.bs_files)['bs']
        # rescale the int values
        bs *= np.log1p(-1.0 / cst.bscl)
        # apply threshold
        bs = np.maximum(np.log(0.01), bs)
        # convert params into udel values
        uvec = np.power(10, cst.params[cst.fixed.bi:cst.fixed.bj])
        # convert udel values into weights relative to initial udel
        bwt = uvec / cst.fixed.u_fix
        # get composite map. rescale the bvals from 0-1000 as ints
        bsx = (1000*np.exp(np.dot(bs, bwt))).astype(int)

        # merge b values that overlap
        bs_un, sg_un = [], []
        last_b = bsx[0]
        last_s = sg[0]
        for (b, s) in izip(bsx[1:], sg[1:]):
            # if the b value hasnt changed, increment the segment
            if b == last_b:
                last_s += s
            # if b switches, save the current segment and bvalue and reset
            else:
                bs_un.append(last_b)
                sg_un.append(last_s)
                last_b = b
                last_s = s

        # add the final segment
        bs_un.append(last_b)
        sg_un.append(last_s)

        # convert to an array
        bs_un = np.array(bs_un)
        sg_un = np.array(sg_un)
        print np.sum(sg_un)
        assert np.sum(sg_un) == cst.chlen

        # save bmap
        f_save = bdir + ch + '.bmap.txt'
        np.savetxt(f_save, np.column_stack((bs_un, sg_un)), fmt='%d %d')


if __name__ == '__main__':
    main()
