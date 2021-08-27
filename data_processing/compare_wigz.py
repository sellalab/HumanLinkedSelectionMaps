__author__ = 'davidmurphy'

import numpy as np
from sys import argv
from classes.runstruct import root_dir, chromosome_length
from classes.wigfixblock import WigfixBlock, wig_iterator, np, zopen


def marked_wigz_array(chrom, fwig):
    """convert a wigz file into an array of values scaled to chromsome"""
    # set mask length, init empty array of integer type for rescaled scores
    mlen = chromosome_length(chrom)
    mask = np.ones(mlen, dtype=int) * -1

    for bi in wig_iterator(wig_file=fwig):
        # make block objects from continuous score blocks in the wig file
        blk = WigfixBlock(data=bi.groups())

        # make sure that blocks obey expected rules
        assert blk.chrom == chrom
        assert blk.step == 1

        # binarize scores into 0s and 1s according to pmin, pmax params
        mask[blk.start:blk.end] = blk.phast_score

    return mask


def compare_files(ch, sp1, sp2):
    """compare conservation scores from two alignments"""
    # set format for files
    # fmt = root_dir + '/data/cons/wigz/{}.{}.wig.gz'
    fmt = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz/{s}/{s}.{c}.wig.gz'

    # get arrays of conservation scores for each of the species
    w1 = marked_wigz_array(ch, fmt.format(s=sp1, c=ch))
    w2 = marked_wigz_array(ch, fmt.format(s=sp2, c=ch))

    s1 = ch + ': {}_covered_{}_missing {}'
    s2 = ch + ': both_missing {}'
    print s1.format(sp2, sp1, np.sum((w1 == -1) & (w2 != -1)))
    print s1.format(sp1, sp2, np.sum((w1 != -1) & (w2 == -1)))
    print s2.format(np.sum((w1==-1) & (w2==-1)))


def measure_intersection(ch):
    """intersect all cons scores and see how many sites are covered"""
    fmt = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz/{s}/{s}.{c}.wig.gz'
    spcs = 'ape primate prosimian euarchontoglires laurasiatheria mammal fish'

    # get fish first (overall most covered bases) and count missing
    w1 = marked_wigz_array(ch, fmt.format(s='fish', c=ch))
    init_miss = np.sum(w1 == -1)

    # loop over remaining species and flip -1s in current to -1 in fish
    for sp in spcs.split():
        w2 = marked_wigz_array(ch, fmt.format(s=sp, c=ch))
        msk = (w2 == -1)
        w1[msk] = -1

    # get final missing count
    final_miss = np.sum(w1 == -1)

    print '{}_init_vs_final_missing {} {}'.format(ch, init_miss, final_miss)


def main():
    if len(argv) != 4:
        print 'usage: compare_wigz <ch> <sp1> <sp2>'
    ch, sp1, sp2 = argv[1:]
    # compare_files(ch, sp1, sp2)
    measure_intersection(ch)

if __name__ == '__main__':
    main()
