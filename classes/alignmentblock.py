import re
import numpy as np
from gzip import open as zopen

__author__ = 'davidmurphy'


def axt_iterator(axt_file):
    """generator that yields (summary, primary, aligned) from axt blocks"""
    # the axt data block pattern
    re_block = re.compile('(\d+ chr\w+ \d+ \d+ chr\w+ \d+ \d+ '
                          '[-+] -?\d+)\n([\w-]+)\n([\w-]+)\n')

    # open the file and read to memory:
    with zopen(axt_file, 'r') as f:
        axt = f.read()

    # block iter returns blocks one at a time from the axt file
    block_iter = re_block.finditer(axt)

    for block in block_iter:
        yield block


class AlignmentBlock(object):
    """
    A class for parsing alignment blocks from UCSC .axt alignment files
    ===================================================================

    Structure
    ---------
    Each alignment block in an axt file contains three lines: a summary
    line and 2 sequence lines. Blocks are separated from one another by
    blank lines.

    1. Summary line
    ---------------

    0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500

    The summary line contains chromosomal position and size information
    about the alignment. It consists of 9 required fields:

    + Alignment number -- The alignment numbering starts with 0 and increments
      by 1, i.e. the first alignment in a file is numbered 0, the next 1, etc.
    + Chromosome (primary organism)
    + Alignment start (primary organism) -- The first base is numbered 1.
    + Alignment end (primary organism) -- The end base is included.
    + Chromosome (aligning organism)
    + Alignment start (aligning organism)
    + Alignment end (aligning organism)
    + Strand (aligning organism) -- If the strand value is "-", the
      values of the aligning organism's start and end fields are relative
      to the reverse-complemented coordinates of its chromosome.
    + Blastz score -- Different blastz scoring matrices are used for
      different organisms. See the README.txt file in the alignments
      directory for scoring information specific to a pair of alignments.

    2. & 3. Sequence lines
    ----------------------

    TCAGCTCATAAATCACCTCCTGCCACAAGCCTGGCCTGGTCCCAGGAGAGTGTCCAGGCTCAGA
    TCTGTTCATAAACCACCTGCCATGACAAGCCTGGCCTGTTCCCAAGACAATGTCCAGGCTCAGA

    The sequence lines contain the sequence of the primary assembly (line 2)
    and aligning assembly (line 3) with inserts. Repeats are indicated by
    lower-case letters.

    (source: http://genome.ucsc.edu/goldenPath/help/axt.html)
    """
    def __init__(self, axt_block):
        """
        Parses the summary line of an axt alignment block and stores
        each alignment string
        """
        # load summary as list
        self._summary = axt_block[0].split()

        # load primary and aligned sequences as uppercase string arrays
        self.primary = np.array(list(axt_block[1].upper()), dtype='S1')
        self.aligned = np.array(list(axt_block[2].upper()), dtype='S1')

        # use primary coordinate system by masking primary gaps in each sequence
        igap = self.primary != '-'
        self.primary, self.aligned = self.primary[igap], self.aligned[igap]

    def __len__(self):
        """the length of the alignment block's reference sequence"""
        return self.ref_end - self.ref_start

    def filtered_submask(self, fpairs):
        """return submask with filter pair sites turned to 0s"""
        # store local copy of submask
        mask = self.submask
        # get indices of mismatch pairs
        sidx = np.where(mask == 2)[0]
        # find mismatch pairs that match filter pairs
        fmsk = np.in1d(self.pairs[sidx], fpairs)
        # flip filter matching pairs to 0s in the submask
        mask[sidx[fmsk]] = 0

        return mask

    @property
    def submask(self):
        """return triplet 012 coded mask where 0=missing; 1=match; 2=mismatch"""
        # initialize all sites as matches (1)
        mask = np.ones(shape=len(self), dtype='u1')

        # recode mismatch sites (2)
        mask[self.primary != self.aligned] = 2

        # recode missing sites if either sequences is missing data (0)
        bases = ['A', 'C', 'G', 'T']
        mask[~np.in1d(self.primary, bases) | ~np.in1d(self.aligned, bases)] = 0

        return mask

    @property
    def pairs(self):
        """return primary/aligned pairs from each sequence"""
        return np.core.defchararray.add(self.primary, self.aligned)

    @property
    def alignment_number(self):
        return int(self._summary[0])

    @property
    def ref_chrom(self):
        return self._summary[1]

    @property
    def ref_start(self):
        return int(self._summary[2]) - 1  # convert to 0-based system

    @property
    def ref_end(self):
        return int(self._summary[3])

    @property
    def align_chrom(self):
        return self._summary[4]

    @property
    def align_start(self):
        return int(self._summary[5])

    @property
    def align_end(self):
        return int(self._summary[6])

    @property
    def align_strand(self):
        return self._summary[7]

    @property
    def blastz_score(self):
        return self._summary[8]
