__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from collections import Counter
from classes.runstruct import ChromStruct, RunStruct, root_dir
from data_processing.data_tools import mask_segments, binary_mask_segments, \
    conservation_mask, dist_percentile
from data_processing.compare_wigz import marked_wigz_array
from classes.geneticmap import GeneticMap
from data_processing.parse_cadd import get_gmap_mask

species = ['ape', 'primate', 'prosimian', 'euarchontoglires', 'laurasiatheria',
           'mammal','fish']


def save_score_dist(score_dist, fout):
    with open(fout, 'w') as f:
        for k in sorted(score_dist.keys()):
            line = '{} {}\n'.format(k, score_dist[k])
            f.write(line)


def clean_scores(ch, neut_ref='euarchontoglires'):
    """
    create a set fo cleaned conservation scores for each alignment such
    that only bases where there is a score from every alignment are given
    scores.

    steps:
    ------
    1. set "neutral" sites to 0
    2. take intersection of sites that have a score for each alignment
    3. save conservation scores for these sites
    4. save distribution of scores from the trimmed set of sites
    """
    # set chromstruct and file path
    cst = ChromStruct(chrom=ch)
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    # set path to genetic map filtered neutmask
    nmsk_pth = root_dir + '/data/nmsk'
    f_nmsk = nmsk_pth + '/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'.format(ch)

    # set input/output file path formats
    fmt_wigz = pth + '/{cn}/{cn}.{c}.wig.gz'
    # fmt_arra = pth + '/{cn}/{cn}.{c}.scores.npz'
    fmt_arra = pth + '/{cn}/{cn}.{c}.gmask.scores.npz'
    # fmt_dist = pth + '/{cn}/{cn}.{c}.clean.score.dist.txt'
    fmt_dist = pth + '/{cn}/{cn}.{c}.gmask.score.dist.txt'

    # load neut mask for default reference alignment=euarchontoglires
    # fneut = fmt_wigz.format(cn=neut_ref, c=ch)
    # fmt = '{p}/{cn}/{cn}.filtered.scores.dist.txt'
    # fdist = fmt.format(p=pth, cn=neut_ref)
    # pmax = dist_percentile(fdist, 0.35)
    # nmsk = conservation_mask(fneut, ch, 0, pmax)

    # load neutmask from file
    nmsk = np.load(f_nmsk)['neutmask']
    # create a mask for sites with scores in all alignments
    # imsk = np.ones(shape=cst.chlen, dtype=bool)
    gmsk = get_gmap_mask(ch)

    # process conservation scores from each alignment
    # scores_dict = {}
    counts_dict = {}
    for cons in species:
        fwigz = fmt_wigz.format(p=pth, cn=cons, c=ch)
        cscores = marked_wigz_array(ch, fwigz)
        # count initial scores
        counts_dict[cons] = [np.sum(cscores != -1)]
        # set neutral positions to zero
        cscores[nmsk] = 0
        # set genetic map masked positions to 0
        cscores[gmsk] = 0
        # count scores after neutral switch
        counts_dict[cons].append(np.sum(cscores != -1))
        # # set sites missing scores (-1s) to False in the score mask
        # imsk &= (cscores != -1)
        # put scores in the dictionary
        # scores_dict[cons] = cscores
        # get score distribution and save
        score_dist = Counter(cscores)
        save_score_dist(score_dist, fmt_dist.format(cn=cons, c=ch))
        # save cleaned conservation scores
        np.savez_compressed(fmt_arra.format(cn=cons, c=ch), p=cscores)

    # # remove sites without scores across alignments using score mask
    # for cons in species:
    #     # scores_dict[cons][~imsk] = -1
    #     # cscores = scores_dict[cons]
    #     # count scores after inersection
    #     # counts_dict[cons].append(np.sum(cscores != -1))
    #     # get score distribution and save
    #     score_dist = Counter(cscores)
    #     save_score_dist(score_dist, fmt_dist.format(cn=cons, c=ch))
    #     # save cleaned conservation scores
    #     np.savez_compressed(fmt_arra.format(cn=cons, c=ch), p=cscores)

    # # print scores before and after filters
    # for cons in species:
    #     msg = '{} {} {} {}'.format(cons, *counts_dict[cons])
    #     print msg

    return None


def clean_conserved_segs(ch, spec, pcons):
    """conserved segments built from cleaned phastcons scores"""
    # set chromstruct and file path
    cst = ChromStruct(chrom=ch)
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'

    # set input/output file path formats and labels
    # anno = '{}_cons{}_clean'.format(spec, pcons)
    # anno = '{}_cons{}_gmask'.format(spec, pcons)
    anno = '{}_cons{}_new'.format(spec, pcons)  # June 2021 update

    bdir = cst.data + '/bsanno/{}'.format(anno)
    # farra = pth + '/{cn}/{cn}.{c}.scores.npz'.format(cn=spec, c=ch)
    farra = pth + '/{cn}/{cn}.{c}.gmask.scores.npz'.format(cn=spec, c=ch)
    # fdist = pth + '/{cn}/{cn}.{c}.clean.score.dist.txt'.format(cn=spec, c=ch)
    fdist = pth + '/{cn}/{cn}.{c}.gmask.score.dist.txt'.format(cn=spec, c=ch)
    fsegs = cst.data + '/bsanno/{an}/{c}.{an}.bed'.format(an=anno, c=ch)

    # create output path if it does not exist
    if not os.path.isdir(bdir):
        os.mkdir(bdir)

    # get the minimum conservation score for the given percentile
    pmin = dist_percentile(fdist, float(pcons) * 0.01)

    # get conservation mask at the given percentage conserved
    cmsk = (np.load(farra)['p'] >= pmin)
    print 'conserved_sites {}'.format(cmsk.sum())

    # convert mask to segments and write to file
    segs = binary_mask_segments(cmsk)
    with open(fsegs, 'w') as f:
        for (start, end) in segs:
            f.write('{} {} {}\n'.format(ch, start, end))

    return None


def clean_conserved_seg_range(ch, spec, pmin, pmax):
    """conserved segments built from cleaned phastcons scores"""
    # set chromstruct and file path
    cst = ChromStruct(chrom=ch)
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'

    # set input/output file path formats and labels
    anno = '{}_cons{}_{}_clean'.format(spec, pmin, pmax)
    bdir = cst.data + '/bsanno/{}'.format(anno)
    farra = pth + '/{cn}/{cn}.{c}.scores.npz'.format(cn=spec, c=ch)
    fdist = pth + '/{cn}/{cn}.{c}.clean.score.dist.txt'.format(cn=spec, c=ch)
    fsegs = cst.data + '/bsanno/{an}/{c}.{an}.bed'.format(an=anno, c=ch)

    # create output path if it does not exist
    if not os.path.isdir(bdir):
        os.mkdir(bdir)

    # get the minimum conservation score for the given percentile
    smin = dist_percentile(fdist, float(pmin) * 0.01)
    smax = dist_percentile(fdist, float(pmax) * 0.01)

    # get conservation mask at the given percentage conserved
    arra = np.load(farra)['p']
    cmsk = (arra >= smin) & (arra <= smax)
    print 'conserved_sites {}'.format(cmsk.sum())

    # convert mask to segments and write to file
    segs = binary_mask_segments(cmsk)
    with open(fsegs, 'w') as f:
        for (start, end) in segs:
            f.write('{} {} {}\n'.format(ch, start, end))

    return None


def main():
    # if len(argv) != 2:
    #     print 'usage: clean_phast_scores <ch>'
    # ch = argv[1]
    # clean_scores(ch)
    if len(argv) == 2:
        ch = argv[1]
        clean_scores(ch)
    elif len(argv) == 4:
        ch, spec, pcons = argv[1:]
        clean_conserved_segs(ch, spec, pcons)
    # elif len(argv) == 5:
    #     ch, spec, pmin, pmax = argv[1:]
    #     clean_conserved_seg_range(ch, spec, pmin, pmax)
    else:
        print 'usage_1: clean_scores <ch>'
        print 'usage_2: clean_conserved_segs <ch> <spec> <pcons>'


if __name__ == '__main__':
    main()
