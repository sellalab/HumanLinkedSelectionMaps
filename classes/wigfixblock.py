from gzip import open as zopen
import numpy as np
import re


__author__ = "davidmurphy"


def wig_iterator(wig_file):
    """a generator that returns WigFixBlocks from a wigFix file"""
    # regex for wigfix block format [chrom, start, stepsize, scores]
    wig = re.compile('fixedStep chrom=(chr[12]?\d|chr[YXM]) start=(\d+) '
                     'step=(\d+)\n((?:-?\d\.\d{3}\n)+)', re.MULTILINE)

    # read in wigFix file as string
    with zopen(wig_file, 'r') as f:
        cons = f.read()

    # block_iter returns contiguous blocks of wig-formatted scores
    block_iter = wig.finditer(cons)

    # iterate over all blocks captured by regex to end of file
    for block in block_iter:
        yield block


class WigfixBlock(object):
    """
    This class represents one block of data in wiggle table fixed-step format
    used by UCSC for phastCons, phyloP scores, etc.

    DOCUMENTATION:
    --------------

    File Format (assemblies released Nov. 2004 and later)

    When uncompressed, the file contains a declaration line and one column of
    data in wiggle table fixed-step format:

      fixedStep chrom=scaffold_1 start=3462 step=1
      0.0978
      0.1588
      0.1919
      0.1948
      0.1684

    1. Declaration line: The declaration line specifies the starting point of
    the data in the assembly. It consists of the following fields:

        fixedStep -- keyword indicating the wiggle track format used to write
                     the data. In fixed step format, the data is single-column
                      with a fixed interval between values.
        chrom -- chromosome or scaffold on which first value is located.
        start -- position of first value on chromosome or scaffold specified
                 by chrom. NOTE: Unlike most Genome Browser coordinates, these
                 are one-based.
        step -- size of the interval (in bases) between values.

    A new declaration line is inserted in the file when the chrom value
    changes, when a gap is encountered (requiring a new start value), or
    when the step interval changes.

    2. Data lines: The first data value below the header shows the score
    corresponding to the position specified in the header. Subsequent score
    values step along the assembly in one-base intervals. The score shows the
    posterior probability that phastCons's phylogenetic hidden Markov model
    (HMM) is in its most-conserved state at that base position.

    (from: http://genome.ucsc.edu/goldenPath/help/phastCons.html)
    """

    def __init__(self, data, pmin=0.0, pmax=1.0):
        """
        Initialize a block from wig fix file that has been broken into
        components with regex
        :param data: A block of data in wig fix format
        :param pmin: minimum score to accept
        :param pmax: max score to accept
        """

        # data should split into 4 regex match groups to use with this class
        assert len(data) == 4

        # initialize member variables from the 4 fields, split scores on '\n'
        self.chrom = data[0]
        self.start = int(data[1])
        self.step = int(data[2])
        self.scores = np.array(data[3].split('\n')[:-1], dtype='f8')
        self.end = self.start + len(self.scores)

        # set internal pmin, pmax for the block
        self.pmin = pmin
        self.pmax = pmax

    @property
    def mask(self):
        """return True if pmin <= score <= pmax else False (boolean)"""
        score_mask = (self.scores >= self.pmin) & (self.scores <= self.pmax)
        return score_mask

    @property
    def phast_score(self):
        """return unit16 scores array on [0, 1000] scale"""
        return (self.scores * 1000).astype('u2')


