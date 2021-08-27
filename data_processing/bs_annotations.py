__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv
from collections import Counter
from classes.runstruct import ChromStruct, RunStruct
from data_processing.data_tools import mask_segments, binary_mask_segments, \
    conservation_mask, dist_percentile, wigz_array, get_filter_mask
from classes.geneticmap import GeneticMap


def get_filtered_cmask(ch, cons, pct, ncons, npct):
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    # path to conserved distribution and wigz file
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)
    # fdist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=cons)
    fdist = '{p}/{cn}/{cn}.filtered.scores.dist.txt'.format(p=pth, cn=cons)
    # path to conserved distribution and neutral wigz file
    fnwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=ncons, c=ch)
    fndist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=ncons)

    # get maximum neutral score value for neutral percentile
    pmax = dist_percentile(fndist, npct)
    # get neutral mask at the cutoff
    nmsk = conservation_mask(fnwigz, ch, 0.0, pmax)

    # get minimum conserved value for the given percentile
    pmin = dist_percentile(fdist, pct)
    # get conserved status mask for values in the selected percentile
    cmsk = conservation_mask(fwigz, ch, pmin, 1.0)

    # take the intersection of conserved and NOT neutral
    filtered_cmsk = (cmsk & ~nmsk).astype('u1')

    return filtered_cmsk


def exnex_conserved_segments(ch, anno):
    """split conserved segments for some pct into exonic, nonexonic"""
    # set output format
    ofmt = '{c} {start} {end}\n'

    # create cst for loading and saving files
    cst = ChromStruct(ch)

    # get annotation mask
    fsegs = cst.bs_target(anno)
    c_seg = np.loadtxt(fsegs, usecols=(1, 2), dtype=int)
    fmask = np.zeros(cst.chlen, dtype=bool)
    fmask = mask_segments(fmask, c_seg, flipoff=False)

    # set the merged coordinate path
    pth = cst.data + '/coords/nr'
    # set exonic file path
    f_x = '{}/exon/{}_knownGene_nonredundant_exonSegments.bed'.format(pth, ch)

    # get a mask with 1s for exonic, 0s for nonexonic
    x_seg = np.loadtxt(f_x, usecols=(1, 2), dtype=int)
    xmask = np.zeros(cst.chlen, dtype=bool)
    xmask = mask_segments(xmask, x_seg, flipoff=False)

    # build conserved exonic filepaths and filenames
    x_anno = anno + '_exonic'
    xpth = '{}/bsanno/{}'.format(cst.data, x_anno)
    # create path to exonic if it doesnt exist
    if not os.path.isdir(xpth):
        os.mkdir(xpth)
    # create save file for exonic
    f_xcon = cst.bs_target(x_anno)
    # create segments of conserved exonic (complement xmask)
    xcon = binary_mask_segments((fmask & xmask))
    with open(f_xcon, 'w') as f:
        for (i, j) in xcon:
            f.write(ofmt.format(c=ch, start=i, end=j))

    # create conserved nonexonic filepaths and filenames
    nx_anno = anno + '_nonexonic'
    nxpth = '{}/bsanno/{}'.format(cst.data, nx_anno)
    # create path to nonexonic if it doesnt exist
    if not os.path.isdir(nxpth):
        os.mkdir(nxpth)
    # create save file for nonexonic
    f_nxcon = cst.bs_target(nx_anno)
    # build segments of conserved nonexonic
    nxcon = binary_mask_segments((fmask & ~xmask))
    with open(f_nxcon, 'w') as f:
        for (i, j) in nxcon:
            f.write(ofmt.format(c=ch, start=i, end=j))

    return None


def cadd_exnex_counts(ch):
    """count total CADD segments and the subset of those that are exonic"""
    # create cst for loading and saving files
    cst = ChromStruct(ch)

    # use CLEAN segments
    anno = 'cadd94'
    fsegs = cst.data + '/bsanno/{an}/{c}.{an}.bed'.format(an=anno, c=ch)

    # get filtered cons mask
    c_seg = np.loadtxt(fsegs, usecols=(1, 2), dtype=int)
    fmask = np.zeros(cst.chlen, dtype=bool)
    fmask = mask_segments(fmask, c_seg, flipoff=False)
    tot_sites = np.sum(fmask)

    # set the merged coordinate path
    pth = cst.data + '/coords/nr'
    # set exonic file path
    f_x = '{}/exon/{}_knownGene_nonredundant_exonSegments.bed'.format(pth, ch)

    # get a mask with 1s for exonic, 0s for nonexonic
    x_seg = np.loadtxt(f_x, usecols=(1, 2), dtype=int)
    xmask = np.zeros(cst.chlen, dtype=bool)
    xmask = mask_segments(xmask, x_seg, flipoff=False)
    tot_exon = np.sum((xmask & fmask))

    # print 'tot/exon {} {}'.format(tot_sites, tot_exon)

    return tot_sites, tot_exon


def main():
    if len(argv) != 3:
        print 'usage: bs_annotations <chrom> <anno>'
        exit(1)

    ch = argv[1]
    anno = argv[2]
    exnex_conserved_segments(ch, anno)


if __name__ == '__main__':
    main()
