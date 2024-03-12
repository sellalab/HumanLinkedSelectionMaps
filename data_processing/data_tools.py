# coding=utf-8
import os
from datetime import timedelta
from shutil import copyfileobj
from collections import defaultdict
from classes.knowngene import complement_strand
from classes.wigfixblock import WigfixBlock, wig_iterator, np, zopen
from classes.alignmentblock import AlignmentBlock, axt_iterator, re

__author__ = 'davidmurphy'


goog_dir = '/Users/davidmurphy/GoogleDrive'
if os.getcwd().startswith('/Users/davidmurphy'):
    root_dir = '{}/linked_selection'.format(goog_dir)
else:
    root_dir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection'


def calculate_error(x_xct, x_est, rel=True, fabs=True):
    """
    Calculate the additive/relative error by comparing set of estimates with
    corresponding exact values for some quantity.

    Additive error = (x_xct - x_est). Relative error = (x_xct - x_est) / x_xct.
    Error terms maybe be asolute value or signed.

    :param x_xct: vector of exact values
    :param x_est: vector of estimates for each quantity in x_xct
    :param rel: use relative error if flagged
    :param fabs: calculate absolute difference |x_xct-x_est| if flagged
    :return:
    """
    error = x_xct - x_est  # get difference
    if rel:
        error /= x_xct  # scale error terms to x_xct
    if fabs:
        error = np.abs(error)  # get absolute error values

    return error


def randint_unique(size, end, start=0):
    """
    Generate a random, unique sample of integers from the range [0, length).

    This function uses numpy.random.randint as a much faster alternative to
    numpy.random.choice (without replacement) for large arrays by creating a
    random sample of integers that can be used as indices to "choose" elements
    of an array instead.

    numpy.random.randint chooses integers with replacement, so the number of
    unique integers in a randint sample may be less than the sample size (i.e.,
    len(current_sample) < size). While (len(current_sample) < size), samples of
    2 * (size - current_size) are concatenated to current_sample and sorted into
    a unique set. This process is repeated until len(current_sample) >= size and
    then the first [0:size] elements of current_sample are returned.

    :param size: sample size
    :param end: max_integer + 1 for sample_rage = [start, end)
    :param start: min_integer for sample_rage = [start, end)
    :return s: sorted sample of unique integers from the line [0, length)
    """
    # convert args to int in case float inputs used
    size, end, start = map(int, (size, end, start))

    # (required to prevent infinitely while loop)
    if not (size <= end - start):
        msg = 'range [{}, {}] too small to accomodate {} random integers'
        raise ValueError(msg.format(start, end, size))

    # generate initial sample of indices
    s = np.unique(np.random.randint(start, end, size))

    # while sample < size: sample additional integers, keep [0:size] uniques
    while s.size < size:
        n = size - s.size
        s_n = np.unique(np.random.randint(start, end, 2*n))
        s = np.unique(np.concatenate((s, s_n)))[0:size]

    # check that sample is unique and sorted
    assert np.all(s[1:] > s[:-1])

    return s


def chromosome_length(chrom):
    """return the integer chromosome length in bp for chrom"""
    chr_length = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
                  'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
                  'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431,
                  'chr10': 135534747, 'chr11': 135006516,
                  'chr12': 133851895, 'chr13': 115169878,
                  'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
                  'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983,
                  'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566,
                  'chrY': 59373566, 'chrX': 155270560,
                  # dmel chroms
                  '2L': 23011544, '2R': 21146708, '3L': 24543557,
                  '3R': 27905053}

    return chr_length[chrom]


def get_filter_mask(chrom, filter_file):
    """create a binary mask from list of bed format filter files"""

    # initialize the mask to "True" at each position
    filter_mask = np.ones(chromosome_length(chrom), dtype=bool)

    # flip positions to "False" at filter file coordinates
    with open(filter_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                fpath = line.strip().format(chrom)
                segments = np.loadtxt(fpath, dtype=int, usecols=(1, 2))
                mask_segments(filter_mask, segments)

    return filter_mask


def binary_mask_segments(mask):
    """
    find [start, end] indices for groups of 1 or more 1s in a binary {01}
    mask
    """
    # 1. create an array from the neutral mask with buffering 0s on either end
    zarr = np.concatenate(([0], mask, [0]))
    # 2. where mask != mask shifted by 1 delineate neutral segments edges
    edges = np.where(zarr[1:] != zarr[:-1])[0]
    # must be an equal number of starts and ends so len(edges) must be even
    assert len(edges) % 2 == 0
    # 3. reshape results (C-style) to get [start, end] coords of each segment
    edges = edges.reshape((len(edges) / 2), 2)

    return edges


def mask_segments(mask, segments, flipoff=True, returnmask=True):
    """
    Return mask with segments flipped to 0 or new segments after removing masked
    regions
    :param mask: chromosome spanning binary mask
    :param segments: 0-based annotation coordinates on the same chromosome
    :param flipoff: flip mask to False in segments if returning mask; map mask
    values onto segment regions if
    returning the segments (default). If flipoff is False, flip mask to True in
    segments if returning mask;
    map mask inverse onto segment regions when returning segments.
    :param returnmask: return mask after segment intersect when true (default),
    otherwise return new segments
    """
    # return the mask after altering regions in segments
    for (i, j) in segments:
        mask[i:j] = 0 if flipoff else 1
    return mask if returnmask else binary_mask_segments(mask)


def sequence_mask(fcall, chrom, call_key='P'):
    """
    return a binary mask from calls, where passing=1, others=0
    :param fcall: path to sequencing call mask file
    :param chrom: name of masked chromosome
    :param call_key: symbol that indicates that a position passed filters
    (default, used by 1000 Genomes, is 'P')
    :return passed: a boolean array where called=True, error=False
    """
    # string -> list -> char array -> boolean mask: P=True, ~P=False
    with zopen(fcall, 'r') as m:
        m.readline()  # skip header line
        chars = np.array(list(m.read().replace('\n', '')))
    assert chars.dtype == '|S1'
    assert chars.size == chromosome_length(chrom)
    passed = (chars == call_key)

    return passed


def fasta_array(chrom, fasta_file):
    """
    load fasta file and convert to numpy string array
    :param chrom: chromosome associated with the fasta file
    :param fasta_file: path to fasta file
    :return farr: fasta data as numpy string array the size of the chrom
    """
    # check for correct chrom in file name
    assert chrom in fasta_file

    # convert fasta to single line upcase string, convert string to array
    with zopen(fasta_file, 'r') as f:
        # skip header line
        f.readline()
        farr = np.array(list(f.read().upper().replace('\n', '')))

    # check for correct size
    assert farr.size == chromosome_length(chrom)

    return farr


def dist_percentile(fdist, pct):
    """return cutoff value for distribution of scores at given percentile"""
    # load distribution of scores from the file
    p, n = np.loadtxt(fdist).T
    # # rescale scores as ints from 0-1000
    # p = (p*1000).astype('u2')
    # rescale counts to cumulative fraction
    n /= n.sum()
    n = np.cumsum(n)
    # get the index where cumulative fraction is >= pct
    i = np.searchsorted(n, pct)
    # adjust for possible overshoot in searchsorted
    i = min(i, len(p)-1)

    # return the p value at the select percentile
    return p[i]


def wigz_array(chrom, fwig, use_neg1=False):
    """convert a wigz file into an array of values scaled to chromsome"""
    # set mask length, init empty array of integer type for rescaled scores
    mlen = chromosome_length(chrom)
    if use_neg1:
        mask = np.ones(mlen, dtype='i2') * -1
    else:
        mask = np.zeros(mlen, dtype='u2')


    for bi in wig_iterator(wig_file=fwig):
        # make block objects from continuous score blocks in the wig file
        blk = WigfixBlock(data=bi.groups())

        # make sure that blocks obey expected rules
        assert blk.chrom == chrom
        assert blk.step == 1

        # binarize scores into 0s and 1s according to pmin, pmax params
        mask[blk.start:blk.end] = blk.phast_score

    return mask


def cpg_mask(chrom, ref):
    """
    Return 0-based positions of reference strand Cs of CpG sites.
    :param chrom: chromosome for CpG mask
    :param ref: reference genome file
    :return: reference strand CpG C positions (0-based)
    """
    # turn fasta file into string array
    farr = fasta_array(chrom, ref)

    # locate all C positions and all G positions
    c = np.where(farr == 'C')[0]
    g = np.where(farr == 'G')[0]

    # locate C positions followed by G positions
    cpg = c[np.in1d(c, g-1)]

    return cpg


def cg_hypermutable():
    """return coordinates for CG hypermutable regions"""
    f_cg = root_dir + '/data/coords/CG_hypermutable.txt'
    # create a dict of Mb blocks with CG hypermutability
    cg_dict = defaultdict(list)
    with open(f_cg, 'r') as f:
        f.next()  # skip header
        for line in f:
            ch, mb, tf = line.strip('\n').split()
            # the T/F indicator is 1 for CG hypermutable segments
            if tf == '1':
                # the Mb given is the upper edge of the segment
                mb_j = int(mb) * 10**6
                mb_i = mb_j - 10**6
                seg = (mb_i, mb_j)
                # add the segment to the list of segments for that chrom
                cg_dict[ch].append(seg)

    return cg_dict


def conservation_mask(fwig, chrom, pmin=0.0, pmax=0.0, mtype=bool, fout=None):
    """
    A function to read phastCons wiggle fixedstep formatted files and
    create a mask based on sites falling within a certain conservation
    score range. The mask is a simple binary mask with 1 for each position
    where the wigfix conservation file score is within [pmin, pmax]. The
    program can be used both with phastCons and phyloP files downloaded
    from UCSC. Missing sites are recorded 0s by default
    :param fwig: the wig file to use for getting cons/ncons regions
    :param chrom: name of masked chromosome
    :param pmin: miniumum acceptable score
    :param pmax: maximum acceptable score
    :param mtype: mask data type (e.g., bool, int, str)
    :param fout: optional file path to store finished mask
    """
    # set mask length
    mlen = chromosome_length(chrom)
    mask = np.zeros(mlen, dtype=np.bool)

    for bi in wig_iterator(wig_file=fwig):
        # make block objects from continuous score blocks in the wig file
        blk = WigfixBlock(data=bi.groups(), pmin=pmin, pmax=pmax)

        # make sure that blocks obey expected rules
        assert blk.chrom == chrom
        assert blk.step == 1

        # binarize scores into 0s and 1s according to pmin, pmax params
        mask[blk.start:blk.end] = blk.mask

    # update the mask to specified data type
    mask = mask.astype(mtype)

    # save to filename, if given
    if fout is not None:
        np.savez_compressed(fout, cmask=mask)

    return mask


def substitution_mask(faxt, chrom, fpairs=None, fout=None):
    """
    Read an alignment file in .axt format, extract each block and then use
    the blocks to construct a mask file with the following score code:
        - 0: either no data on one or both sequences, or to pad a spatial
             gap between successive blocks
        - 1: true ACGT base in both sequences and seq1 == seq2
        - 2: true ACGT base in both sequences but seq1 != seq2

    :param faxt: optional faxt to process (default to stored
    filename in parent class)
    :param chrom: name of masked chromosome
    :param fpairs: optional set of mismatch pairs to filter out from alignment
    :param fout: a flag that indicates that the substitution mask
    should be saved to the default file location specified in RunStruct
    :return mask: numpy array the length of chrom using the 012 encoding
    described above or None if fout=True
    """
    # set mask length
    mlen = chromosome_length(chrom)
    mask = np.zeros(mlen, dtype='u1')

    # iterate through the alignment and build the mask
    num = -1  # block numbers given in the file
    end = 0  # current end pos
    for bi in axt_iterator(axt_file=faxt):
        # take the match object and make it into an AlignmentBlock object
        blk = AlignmentBlock(axt_block=bi.groups())

        # idx, pos increase correctly, seq sizes match, ch is correct
        assert blk.ref_chrom == chrom
        assert end <= blk.ref_start
        assert blk.alignment_number == num + 1
        assert len(blk.primary) == len(blk.aligned)

        # add block mask sequence at the block's chromosomal coordinates
        if fpairs is None:
            mask[blk.ref_start:blk.ref_end] = blk.submask
        else:
            # apply optional pair filtering (e.g., biased gene conversion sites)
            mask[blk.ref_start:blk.ref_end] = blk.filtered_submask(fpairs)

        # reset positional indicators
        end = blk.ref_end
        num = blk.alignment_number

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, submask=mask)

    return mask


def substitution_typemask(chrom, faxt, fanc, fout=None):
    # TODO: finish this when the status of X-Aut divergence becomes clear
    """
    For specific mutational classes. Let’s start with:
    (1) CpG>TpG
    (2) C>G
    (3) all types EXCEPT 1&2 above
    (4) C>T

    Include the same pattern on the opposite strand (i.e. C>G includes G>C,
    as is always done, etc). Full list: C>G, C>T, C>A, A>G, A>C, A>T & CpG’s.

    Use the following numeric encoding to stratify matches and mismatches by
    ancestral state:
    0: missing data
    1, 2: are used when at least the primary genome matches the ancestor
    for A/T, C/G, respectively (only count subs on the primary lineage)
    3, 4, 5: A>C/T>G, A>G/T>C, A>T/T>A
    6, 7, 8: C>A/G>T, C>G/G>C, (NON-CpG) C>T/G>A
    9: CpG>TpG/GpC>ApC


    :param chrom: chromosome used for alignments
    :param faxt: alignments file in .axt format
    :param fanc: ancestor file in .fa format (should correspond to alignment)
    :param fout: optional save path for final mask file
    :return mask: mask encoding each type of substitution detected
    """
    # load up ancestral sequence for the chromosome
    aa = fasta_array(chrom, fanc)

    # set mask length
    mlen = chromosome_length(chrom)
    mask = np.zeros(mlen, dtype='u1')

    # iterate through the alignment and build the mask
    num = -1  # block numbers given in the file
    end = 0  # current end pos
    for bi in axt_iterator(axt_file=faxt):
        # take the match object and make it into an AlignmentBlock object
        blk = AlignmentBlock(axt_block=bi.groups())

        # idx, pos increase correctly, seq sizes match, ch is correct
        assert blk.ref_chrom == chrom
        assert end <= blk.ref_start
        assert blk.alignment_number == num + 1
        assert len(blk.primary) == len(blk.aligned)

        # get the submask for the block to identify passing sites
        submsk = blk.submask()

        # add block mask sequence at the block's chromosomal coordinates
        mask[blk.ref_start:blk.ref_end] = blk.submask()

        # reset positional indicators
        end = blk.ref_end
        num = blk.alignment_number

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, submask=mask)

    return mask


def substitution_counts(neutmask, window, chrom, fout=None):
    """
    Count match/mismatch neutral bases per specifed window size
    :param neutmask: FASTA-style mask encoded by 0s, 1s and 2s with the
    translated as follows:
        0=non-neutral/missing
        1=neutral & matches outgroup
        2=neutral & diverged from outgroup
    neutmask can be created with neutrality_mask function
    :param window: set different window size for counting neutral
    match/mismatch positions
    :param chrom: name of masked chromosome
    :param fout: optional file path to store finished array
    should be saved to the default file  location specified in RunStruct
    :return [window, bases, subs] array for the chrom
    If write flag set to True, files are written with the following format:
        # POS   BASES   SUBS
        100000  5432    834
    POS: upper bound of bin, i.e. rates from 0-1000 will have POS = 1000
    BASES: total number of "1" and "2" in the window, i.e. all valid bases
    SUBS: total number of "2" in the window, just the substitutions
    """
    chlen = chromosome_length(chrom)
    # calculate the total number of windows that will cover the chromosome
    tot_bins = int(np.ceil(1.0 * chlen / window))
    # empty array for [window_end, tot_bases, tot_subs] for each window
    counts_array = np.zeros(shape=(tot_bins, 3), dtype='u4')

    for idx in range(tot_bins):
        i, j = int(idx*window), int((idx + 1)*window)
        # count match (1+3) and mismatch (2+6) in neutmask block
        nblk = neutmask[i:j]
        match = np.sum((nblk == 1) | (nblk == 3))
        mismatch = np.sum((nblk == 2) | (nblk == 6))
        # match, mismatch = np.bincount(neutmask[i:j], minlength=3)[1:]
        # record [upper_edge, tot_matches, tot_mismatches] per window
        counts_array[idx] = min(j, chlen), match + mismatch, mismatch

    # substitution rate estimates should tile the full chromosome
    assert counts_array[-1, 0] == chlen

    if fout:
        np.savez_compressed(fout, subcount=counts_array)

    return counts_array


def neutrality_mask(nconsmask, callmask, submask, isnp,
                    filtmask=None, fout=None):
    """
    A mask of the final set of fully filtered neutral positions.

                      0=non-neutral/missing
                      1=neutral, monomorphic, matches outgroup
                      2=neutral, monomorphic, diverged from outgroup
                      3=neutral, polymorphic, matches outgroup
                      6=neutral, polymorphic, diverged from outgroup

    :param nconsmask: mask of neutral, nonconserved sites -> can be created
                      using conservation mask function
    :type nconsmask: np.ndarray
    :param callmask: mask of polymorphism calls
    :type callmask: np.ndarray
    :param submask: mask of match/mismatch with an outgroup
    :type submask: np.ndarray
    :param isnp: (1-based) SNP positions along the chromosome
    :type isnp: np.ndarray
    :param filtmask: (optional) mask of 1 or more combined filters
    :type filtmask: np.ndarray
    :param fout: optional file path to store finished mask
    :type fout: str
    :return neutmask: FASTA-style mask encoded by 0s, 1s and 2s with the
                      translated as follows:
    :rtype np.ndarray
    """
    # initialize each position to True (i.e. neutral)
    neutmask = np.ones(shape=nconsmask.size, dtype=bool)

    # default mask: 100-vertebrates phastCons score == 0
    neutmask &= nconsmask

    # remove non-passing sites in polymorphism call mask (bool)
    neutmask &= callmask

    # apply optional filtering mask
    if filtmask is not None:
        neutmask &= filtmask

    # encode substitution status at neutral sites (removes nonaligned sites)
    neutmask = neutmask.astype('u1') * submask  # uint8 preserves 012 encoding

    # encode polymorphism status at neutral sites (poly+match=3, poly+div=6)
    neutmask[isnp-1] *= 3

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, neutmask=neutmask)

    return neutmask


def neutral_snps(neutmask, snpcnt, fout=None):
    """
    Create array of (0-based) [positions, sample_size, alt_count] for SNP data.

    :param neutmask: FASTA-style mask encoded by 0s, 1s and 2s with the
    translated as follows:
        0=non-neutral/missing
        1=neutral & matches outgroup
        2=neutral & diverged from outgroup
    neutmask can be created with neutrality_mask function
    :param snpcnt: array of polymorphic sites formatted:
        SNP_POS REF_COUNT ALT_COUNT
    :param fout: optional file path to store finished mask
    :return snps:
    """
    # convert SNP positions to 0-based coordinate system
    pos = snpcnt[:, 0]
    pos -= 1  # NOTE: change numpy.view pos changes snpcnt as well

    # get boolean mask of neutral SNP sites
    isneut = neutmask[pos] > 0
    neutsnps = snpcnt[isneut]

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, neutsnps=neutsnps)

    return neutsnps


def snpcount(fsnp, returnbases=False):
    """
    Get polymorphism data from vcftools formatted .frq.count files
    :param fsnp: snp_count file formatted (after vcftools):
        CHROM POS N_ALLELES N_CHR REF:COUNT ALT:COUNT
    :param returnbases: optionally return the REF/ALT bases
    :return snp_data: [pos, refcount, altcount] or
    [pos, ref, refcount, alt, altcount] if returnbases=True
    ref, refcount, alt, altcount]
    :rtype: tuple
    """
    # match [snp_pos, ref_base, ref_count, alt_base, alt_count]
    if returnbases:
        re_snps = re.compile(
            '[\dX]{1,2}\s(\d{1,9})\s\d+\s\d+\s([ACGT]):(\d+)\s([ACGT]):(\d+)')
        with zopen(fsnp, 'r') as f:
            ar = np.array(re_snps.findall(f.read())).T
        snppos = ar[0].astype('u4')
        rbase = ar[1]
        rcount = ar[2].astype('u4')
        abase = ar[3]
        acount = ar[4].astype('u4')

        return snppos, rbase, rcount, abase, acount

    # match [snp_pos, ref_count, alt_count]
    else:
        re_snps = re.compile(
            '[\dX]{1,2}\s(\d{1,9})\s\d+\s\d+\s[ACGT]:(\d+)\s[ACGT]:(\d+)')
        with zopen(fsnp, 'r') as f:
            ar = np.array(re_snps.findall(f.read().decode('utf-8')), dtype='u4').T
        snppos = ar[0]
        rcount = ar[1]
        acount = ar[2]

        return snppos, rcount, acount


def snp_types(chrom, fsnp, fanc, fout=None):
    """
    Identify mutation types for human SNP based human:chimp ancestral sequence.
    :param chrom: chromosome where SNPs reside
    :param fsnp: SNP counts file
    :param fanc: ancestral reference FASTA file
    :param fout: optional filename for [POS, TYPE] formatted file
    :return snptype: [POS, TYPE] formatted array, where pos is 1-based
    :rtype: np.ndarray
    """
    # check that the chrom is consistent across files
    assert (chrom in fsnp) and (chrom in fanc)

    # get SNP data including ref/alt bases
    snppos, rbase, rcount, abase, acount = snpcount(fsnp, returnbases=True)

    # get ancestral sequence and bases at SNP sites
    aa_seq = fasta_array(chrom, fanc)
    aa = aa_seq[snppos-1]

    # valid sites to analyze, where ancestor matches either ref or alt
    # -> the XOR statement excludes sites where neither ref or alt matches AA
    msk = (aa != rbase) ^ (aa != abase)

    # filter non-valid sites
    snppos, rbase, abase, aa = [ar[msk] for ar in [snppos, rbase, abase, aa]]

    # get sites where ref or alt match the ancestral state
    ridx = (aa == rbase)
    aidx = (aa == abase)

    # create arrays of 2-character anc:der pairs where alt or ref are derived
    asnp = aa[ridx].astype(np.chararray) + abase[ridx].astype(np.chararray)
    rsnp = aa[aidx].astype(np.chararray) + rbase[aidx].astype(np.chararray)

    # combine anc-der pairs into a single array maintaining positional order
    pairs = np.chararray(len(snppos), 2)
    pairs[ridx] = asnp
    pairs[aidx] = rsnp

    # annotate each time of anc-der pair
    snptype = np.zeros(shape=snppos.size, dtype='u1')
    for (i, p) in enumerate('AC AG AT CA CG CT'.split(), 1):
        # find sites for each pair and its reverse complement
        imatch = (pairs == p) | (pairs == complement_strand(p))

        # encode the SNP type at each site
        snptype[imatch] = i

    # recode CT/GA sites with CpG context in the ancestral sequence
    ct_idx = np.where(pairs == 'CT')[0]
    ga_idx = np.where(pairs == 'GA')[0]

    # encode CpG/non-CpG transitions
    snptype[ct_idx[aa_seq[snppos[ct_idx]] == 'G']] = 7
    snptype[ga_idx[aa_seq[snppos[ga_idx]] == 'C']] = 7

    if fout:
        pass

    return np.column_stack((snppos, snptype))


def time2str(time_obj):
    # use a strict hours:minutes:seconds format to store time durations
    assert isinstance(time_obj, timedelta)
    tot_sec = int(time_obj.total_seconds())
    seconds = tot_sec % 60
    minutes = (tot_sec % 3600) / 60
    hours = tot_sec / 3600
    return '{:02}:{:02}:{:02}'.format(hours, minutes, seconds)


def str2time(time_str):
    retime = re.compile('((?P<days>\d+)\sdays*,\s)*(?P<hours>\d{1,2}):'
                        '(?P<minutes>\d{1,2}):(?P<seconds>[\d.]+)')

    # convert "n day(s), hours:minutes:seconds" to a datetime.timedelta object
    if not isinstance(time_str, str):
        msg = 'ERROR: time_str is {}, not <str>'.format(type(time_str))
        raise TypeError(msg)
        # assert isinstance(time_str, str)
    try:
        timedict = retime.search(time_str).groupdict()
        for (k, v) in timedict.items():
            if v:
                timedict[k] = float(v)
            else:
                timedict[k] = 0
        return timedelta(**timedict)
    except ValueError:
        msg = 'Elapsed time format <{}> not recognized'
        raise ValueError(msg.format(time_str))
    except AttributeError:
        msg = 'Elapsed time format <{}> not recognized'
        raise AttributeError(msg.format(time_str))


def f_zip(f, fz=None, rm=False):
    """make a zipped file f (to gzip 'best' equivalent)"""
    # default fz is at the same address with ".gz" suffix added
    if fz is None:
        fz = f + '.gz'
    # create the new compressed file
    with open(f, 'r') as f_in, zopen(fz, 'w') as f_out:
        copyfileobj(f_in, f_out)
    # optional remove original
    if rm:
        os.remove(f)


def f_unzip(fz, f=None, rm=False):
    """make an unzipped file f"""
    # default f is at the same address with ".gz" suffix removed
    if f is None:
        f = fz[:-3]
    # create the decompressed file
    with zopen(fz, 'r') as f_in, open(f, 'w') as f_out:
        copyfileobj(f_in, f_out)
    # optional remove original
    if rm:
        os.remove(fz)


def efficient_datatype(max_value, data_type):
    """return smallest the data type to accomodate the set of values"""
    # number of bits = base_2 exponents of each max value per dtype
    nbits = np.logspace(3, 6, 4, base=2).astype('uint8')

    # signed ints
    if data_type == 'int':
        # if value is negative, use equivalent positive value in terms of bits
        if max_value < 0:
            max_value = abs(max_value) - 1
        # use bins to find the min nbits that can carry max_val
        bins = [2**x - 1 for x in nbits - 1]
        idx = np.searchsorted(bins, max_value)
        # return dtype name if max_value is within range
        if 0 <= idx <= 3:
            return 'int{}'.format(nbits[idx])
        else:
            msg = 'Value <{}> too large for dtype <int>'
            raise ValueError(msg.format(max_value, data_type))

    # analagous processing for unsigned ints
    elif data_type == 'uint':
        if max_value < 0:
            msg = 'Negative values not compatible with dtype <uint>'
            raise ValueError(msg.format(max_value))
        bins = [2**x for x in nbits]
        idx = np.searchsorted(bins, max_value)
        if 0 <= idx <= 3:
            return 'uint{}'.format(nbits[idx])
        else:
            msg = 'Value <{}> too large for dtype <uint>'
            raise ValueError(msg.format(max_value, data_type))

    # boolean
    elif data_type == 'bool':
        return 'bool'

    # floats
    elif data_type == 'float':
        # TODO: handling floats efficiently (for now use float64 for safety)
        return 'float64'

    # other?
    else:
        raise ValueError('Data type <{}> was not recognized'.format(data_type))

