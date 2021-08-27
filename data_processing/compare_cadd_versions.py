__author__ = 'davidmurphy'


import numpy as np
from gzip import open as gzopen
from sys import stderr, stdout, argv
from data_processing.parse_cadd import CADDLine
from classes.runstruct import chromosome_length


cadd_root = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/cadd'
cadd_dir1 = cadd_root + '/cadd_v1.4'
cadd_dir2 = cadd_root + '/cadd_v1.6_without_bstat'
cadd_dir3 = cadd_root + '/cadd_v1.6'


def compare_cadd_files(chrom):
    """compare the old and new cadd files for different top score cutoffs"""
    # get the old CADD scores
    f1 = cadd_dir3 + '/scores/{}.cadd.scores.npz'.format(chrom)
    a1 = np.load(f1)['cadd']

    # get the new CADD scores
    f2 = cadd_dir2 + '/scores/{}.cadd.scores.npz'.format(chrom)
    a2 = np.load(f2)['cadd']

    # for percentiles 2-9, compare top scores
    for pct in range(2, 10):
        # get the cutoff for old scores and mask of sites >= cutoff
        cut1 = np.percentile(a1, 100-pct)
        msk1 = (a1 >= cut1)
        l1 = np.sum(msk1)

        # repeat for new scores
        cut2 = np.percentile(a2, 100-pct)
        msk2 = (a2 >= cut2)
        l2 = np.sum(msk2)

        # count the sites where both masks match
        match12 = np.sum(msk1 & msk2)

        # get the percentage matching
        p1 = 100.0 * match12 / l1
        p2 = 100.0 * match12 / l2

        # print a series of messages
        msg =  'Comparison for {} at {}%\n'.format(chrom, pct)
        msg += 'total sites: v1.6={} v1.6_without_B={}\n'.format(l1, l2)
        msg += 'total matches: v1.6={:.6f} v1.6_without_B={:.6f}\n'.format(p1, p2)
        msg += '----\n'
        stderr.write(msg)
        stdout.flush()


def main():
    if len(argv) != 2:
        print 'usage: compare_cadd_versions <chrom>'
        exit(1)

    chrom = argv[1]
    compare_cadd_files(chrom)


if __name__ == '__main__':
    main()
