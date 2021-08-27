from data_tools import chromosome_length
from classes.runstruct import root_dir
from merge import merge_weighted
from sys import argv
import numpy as np


__author__ = 'davidmurphy'


def fix_gaps(chrom, newlen, bedfile):
    """
    Given Fill out the gaps in lifted over .bed file so that the segments
    sum to the new build chromosome length
    :param chrom: chromosome label
    :param newlen: length of chrom in new build
    :param bedfile: lifted over bed file
    """
    # load the coordinate data and b vals
    coords = np.loadtxt(bedfile, usecols=(1, 2, 3), dtype='int32')
    start, end, bval = coords.T

    # get sorting index for start values and sort
    sidx = np.argsort(start)
    start, end, bval = coords[sidx].T

    # merge sorted coords
    merger = merge_weighted(np.column_stack((start, end)), bval)
    start, end, bval = np.array(list(merger)).T

    # assert that there are no rearrangements
    assert np.all(end > start)
    assert np.all(start[1:] > end[:-1])

    # adjust the first start coordinate to 1
    start[0] = 1

    # adjust the remaining start coordinates to previous end + 1
    start[1:] = end[:-1]+1

    # adjust the final end to the new length
    end[-1] = newlen

    # get the new segment lengths for .bkgd formatting
    segments = end - start + 1

    # assert that the new sum of segments is correct
    assert np.sum(segments) == newlen

    # save new .bkgd format files
    new_bkgdfile = bedfile.replace('.bed', '')
    np.savetxt(new_bkgdfile, np.column_stack((bval, segments)), fmt='%d %d')

    # rejoin the fixed data and convert to str
    data = np.column_stack((start, end, bval)).astype(str)

    # create column of chrom labels for bed formatted file
    stype = 'S{}'.format(len(chrom))
    clabels = np.full(len(data), chrom, dtype=stype)
    data = np.column_stack((clabels, data))

    # re-write the gap-filled file
    new_bedfile = bedfile.replace('.bed', '.gapfilled.bed')
    np.savetxt(new_bedfile, data, fmt='%s %s %s %s')


def main():
    if root_dir.startswith('/Users/davidmurphy/'):
        for chrom in ['chr{}'.format(c) for c in range(3, 23) + ['X']]:
            # chrom = 'chr2'
            bdir = '{}/mcvicker/original_hg18'.format(root_dir)
            bedfile = '{}/{}.hg19.bkgd.bed'.format(bdir, chrom)
            newlen = chromosome_length(chrom)
            fix_gaps(chrom, newlen, bedfile)
            print '{} done.'.format(chrom)

    # else:
    #     chrom = argv[1]
    #     bedfile = argv[2]
    #
    # newlen = chromosome_length(chrom)
    # fix_gaps(chrom, newlen, bedfile)


if __name__ == '__main__':
    main()
