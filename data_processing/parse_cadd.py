__author__ = 'davidmurphy'


import os
import numpy as np
from gzip import open as gzopen
from collections import Counter
from sys import stdin, stderr, stdout, argv
from classes.runstruct import chromosome_length, human_autosomes, ChromStruct
from data_processing.data_tools import binary_mask_segments, dist_percentile
from classes.geneticmap import GeneticMap



# cadd_dir = cadd_dir3
anno_dir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/bsanno'
root_dir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection'
nmsk_pth = root_dir + '/data/nmsk'


# NOT USED!
# def get_neutral_gmask(ch):
#     """get the neutral mask with gmap masking"""
#     # fstr = '/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'.format(ch)
#     fstr = '/{}.euarchontoglires.0.35.gmsk0.0.nmsk.npz'.format(ch)
#     f_nmsk = nmsk_pth + fstr
#     nmsk = np.load(f_nmsk)['neutmask']
#
#     return nmsk


def get_gmap_mask(ch):
    cst = ChromStruct(ch)
    # remove genetic map edges
    gmsk = np.zeros(shape=cst.chlen, dtype=bool)
    # load the gmap
    gmp = GeneticMap(cst.chrom, cst.gmap_files)
    # get genetic map position to mask up to from 5' direction
    gpos5 = gmp.interp_gpos(gmp.viable_5prime) + cst.gmsk
    # get physical position corresponding to the genetic map masking
    pos5 = gmp.interp_pos(gpos5)
    # get the index to maximum non-zero position, set region to False
    gmsk[:pos5] = True
    # get genetic map position to mask starting in the 3' direction
    gpos3 = gmp.interp_gpos(gmp.viable_3prime) - cst.gmsk
    # get physical position corresponding to the genetic map masking
    pos3 = gmp.interp_pos(gpos3)
    # get the index to maximum non-zero position, set region to False
    gmsk[pos3:] = True

    return gmsk


class CADDLine:
    """a structure that holds subset of data from CADD file lines"""
    def __init__(self, input_line):
        """Split line of CADD input into fields"""
        # 1	10001	T	A	0.118631	4.575
        data = input_line.split()
        self.chrom = 'chr{}'.format(data[0])
        self.pos = int(data[1])
        self.phred = float(data[-1])


class CADDReader:
    """
    a class for reading and processing CADD data from the file:
    https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/
    whole_genome_SNVs.tsv.gz
    (more info at: https://cadd.gs.washington.edu/download)
    """
    def __init__(self, cadd_file, cadd_dir):
        self.line_count = 0
        self.pct = 0
        self.scores = Counter()
        self.chrom = 'chr1'
        self.pos = 0
        self.trip = []
        self.arr = self.newarray
        self.cadd_file = cadd_file
        self.cadd_dir = cadd_dir

    def _checkcount(self):
        """routine for printing a progress bar in the log file"""
        # count line and print '.' every 85M lines (approx 1%)
        self.line_count += 1
        if self.line_count % 85000000 == 0:
            self.pct += 1
            if self.pct % 10 == 0:
                msg = '{}%'.format(self.pct)
            else:
                msg = '.'
            stderr.write(msg)
            stdout.flush()

    def _add2dict(self):
        """routine for adding value to the score dictionary"""
        # record max phred score at position to counter
        if len(self.trip) == 3:
            self.scores[max(self.trip)] += 1
            self.trip = []
        else:
            assert len(self.trip) == 0

    def _add2array(self):
        """routine for adding value to current position in array"""
        if len(self.trip) == 3:
            self.arr[self.pos - 1] = max(self.trip)
            self.trip = []
        else:
            assert len(self.trip) == 0

    def _nextchrom(self, caddline):
        """routine for finishing current chrom and starting new one"""
        # check that conditions are corrext
        assert caddline.pos != self.pos
        assert len(self.trip) == 3

        # fill last position in current array
        self.arr[self.pos-1] = max(self.trip)

        # save current array
        fmt_arr = self.cadd_dir + '/scores/{}.cadd.scores.npz'.format(self.chrom)
        np.savez_compressed(fmt_arr, cadd=self.arr)

        # print message
        stderr.write('{} complete.\n'.format(self.chrom))
        stdout.flush()

        # initialize new array, reset chrom, pos, trip
        self.chrom = caddline.chrom
        self.arr = self.newarray
        self.pos = caddline.pos
        self.trip = []

    def _savedict(self):
        # after completing the score count, record districution to file
        fout = self.cadd_dir + '/cadd.scores.txt'
        with open(fout, 'w') as f:
            for score in sorted(self.scores.keys()):
                line = '{} {}\n'.format(score, self.scores[score])
                f.write(line)

    def fetch_data(self):
        f = self.open_file
        for line in f:
            # skip comments
            if line.startswith('#'):
                continue

            # use raw data to initialize CADDLine class
            caddline = CADDLine(line)

            # finish current chrom and start new one
            if caddline.chrom != self.chrom:
                self._nextchrom(caddline)

            if self.pos != caddline.pos:
                self._add2array()

            self.pos = caddline.pos
            self.trip.append(caddline.phred)

        f.close()

    def profile_data(self):
        f = self.open_file
        for line in f:
            # skip comments
            if line.startswith('#'):
                continue

            # monitor progress
            self._checkcount()

            # use raw data to initialize CADDLine class
            caddline = CADDLine(line)

            # add max score at current position to dict, move to next position
            if self.pos != caddline.pos:
                self._add2dict()

            self.pos = caddline.pos
            self.trip.append(caddline.phred)

        self._savedict()


    @property
    def newarray(self):
        return np.zeros(shape=chromosome_length(self.chrom))

    @property
    def open_file(self):
        return gzopen(self.cadd_file, 'r')


class CADDScore:
    """a class to process CADD scores in various ways"""
    def __init__(self, chrom, cadd_dir):
        f_in = cadd_dir + '/scores/{}.cadd.scores.npz'.format(chrom)
        self.scores = np.load(f_in)['cadd']
        self.chrom = chrom
        self.cadd_dir = cadd_dir

    def get_mask(self, pct):
        """return a binary mask for scores over/under threshold"""
        # TODO: don't leave this on the old distribution forever
        # f_dist = cadd_dir + '/cadd.scores.txt'
        # f_dist = cadd_dir + '/cadd.scores.corrected.txt'
        f_dist = self.cadd_dir + '/cadd.scores.gmask.corrected.txt'
        threshold = dist_percentile(f_dist, 0.01*pct)
        # load neutral mask to genetic map masking
        gmsk = get_gmap_mask(self.chrom)
        # set genetic map masked sites to 0 to shift distribution
        self.scores[gmsk] = 0

        return self.scores >= threshold

    def get_segs(self, pct):
        """get segments from scores passing a certain threshold"""
        # create a binary mask based on over/under threshold
        mask = self.get_mask(pct)
        print('{} {}'.format(self.chrom, mask.sum()))
        segs = binary_mask_segments(mask)
        return segs

    def save_segs(self, pct, label):
        """save segments passing threshold to BED file"""
        # make save directory if it does not exist
        f_path = anno_dir + '/cadd{}_gmask_{}'.format(int(pct), label)
        if not os.path.isdir(f_path):
            os.mkdir(f_path)
        # create file save path
        ffmt = f_path + '/{}.cadd{}_gmask_{}.bed'
        f_save = ffmt.format(self.chrom, int(pct), label)
        # get segments and save to file
        segs = self.get_segs(pct)
        with open(f_save, 'w') as f:
            for (start, end) in segs:
                f.write('{} {} {}\n'.format(self.chrom, start, end))

    @property
    def score_dist(self):
        return Counter(self.scores)


def profile_autosomes(cadd_dir):
    """get distribution of autosomal CADD scores"""
    # counter for each CADD score
    score_dict = Counter()

    # get scores from each autosome and add occurrences to counter
    for ch in human_autosomes:
        cadd = CADDScore(ch, cadd_dir)
        # mask out genetic map edge scores
        gmsk = get_gmap_mask(ch)
        cadd.scores[gmsk] = 0
        score_dict += cadd.score_dist
        msg = '{} done.\n'.format(ch)
        stderr.write(msg)
        stdout.flush()

    # print the total number of values seen
    msg = 'total scores: {}\n'.format(sum(score_dict.values()))
    stderr.write(msg)
    stdout.flush()

    # after completing the score count, record distribution to file
    # fout = cadd_dir + '/cadd.scores.corrected.txt'
    fout = cadd_dir + '/cadd.scores.gmask.corrected.txt'
    with open(fout, 'w') as f:
        for score in sorted(score_dict.keys()):
            line = '{} {}\n'.format(score, score_dict[score])
            f.write(line)


def main():
    cadd_root = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/cadd'
    cadd_dir1 = cadd_root + '/cadd_v1.4'  # original
    cadd_dir2 = cadd_root + '/cadd_v1.6_without_bstat'
    cadd_dir3 = cadd_root + '/cadd_v1.6'

    # if len(argv) != 3:
    #     # print 'usage: parse_cadd <cadd_file> <flag>'
    #     print 'usage: parse_cadd <chrom> <pct>'
    #     exit(1)

    # f_in = argv[1]
    # reader = CADDReader(f_in)
    # if argv[2] == 'profile':
    #     reader.profile_data()
    # elif argv[2] == 'fetch':
    #     reader.fetch_data()
    # else:
    #     print 'option <{}> does not exist'.format(argv[2])
    chrom = argv[1]
    pct = float(argv[2])

    # make segments for both sets of scores in one call
    labels = ['v1.6_without_bstat', 'v1.6']
    caddirs = [cadd_dir2, cadd_dir3]
    # for (lbl, cdir) in zip(labels, caddirs):
    CADDScore(chrom, cadd_dir2).save_segs(pct, labels[0])

    # _, _ = argv[1:3]
    # profile_autosomes(cadd_dir)


if __name__ == '__main__':
    main()

