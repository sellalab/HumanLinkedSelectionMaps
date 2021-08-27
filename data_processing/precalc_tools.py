import numpy as np
from itertools import izip
from data_tools import neutral_snps
from data_processing.data_tools import chromosome_length
from functions import count_all_pairs

__author__ = 'davidmurphy'


def open_filetype(fname):
    """a file reader that uses different protocols by filetype"""
    # basic npy file
    if fname.endswith('.npy'):
        data = np.load(fname)

    # npz file with single archive
    elif fname.endswith('.npz'):
        z = np.load(fname)
        k = z.keys()
        assert len(k) == 1

        data = z[k[0]]

    # everything else will be opened as a text file for now
    else:
        data = np.loadtxt(fname)

    return data


def collect_bmaps(chrom, bmap_files):
    """
    collect bmap segments and bmap values into two separate lists for a set
    of bmap files
    :param bmap_files: a list of bmap files
    :param chrom: chromosome name shared by the files
    :return segments, values: two lists of numpy arrays containing segment
    lengths and b values, respectively
    """
    # check that the files are for the same chromosome
    assert all(chrom in f for f in bmap_files)
    # TODO: integrate BkgdMap and BkgdMapReader classes here
    segments, values = [], []
    for f in bmap_files:
        # split B values and segment lengths into separate arrays
        bvals, segs = np.loadtxt(f).T
        # check that bmaps tile the entire chromosome
        assert segs.sum() == chromosome_length(chrom)
        segments.append(segs.astype('u8'))
        values.append(bvals.astype('f8'))

    return segments, values


def collect_lsmaps(chrom, map_files):
    """
    collect BS/CS segments and values into two separate lists
    :param map_files: a list of map files
    :param chrom: chromosome name shared by the files
    :return segments, values: two lists of numpy arrays containing segment
    lengths and cs/bs values, respectively
    """
    # check that the files are for the same chromosome
    assert all(chrom in f for f in map_files)
    segments, values = [], []
    for f in map_files:
        # split values and segment lengths into separate arrays
        vals, segs = open_filetype(f).T
        # check that bmaps tile the entire chromosome
        assert segs.sum() == chromosome_length(chrom)
        segments.append(segs.astype('u8'))
        values.append(vals.astype('f8'))

    return segments, values


def prepare_umap(chrom, sub_counts):
    """mini-function to prepare the umap segments and values"""
    bins, bases, mismatch = sub_counts.astype('f8').T
    # bins should tile total chrom length
    assert bins[-1] == chromosome_length(chrom)

    # get the segment length of each bin
    segs = bins - np.concatenate(([0], bins[:-1]))

    # mask sites with too little data
    mask = (bases >= 100) & (mismatch > 0)

    # estimate rates in the remaining windows
    rates = np.zeros(shape=len(bins))
    rates[mask] = mismatch[mask] / bases[mask]

    # interpolate rate estimates for masked bins from neighboring bins
    rates[~mask] = np.interp(bins[~mask], bins[mask], rates[mask])

    return segs.astype('u8'), rates


def prepare_umap_2(chrom, subrate):
    """mini-function to prepare the umap segments and values"""
    bins, srate = subrate.astype('f8').T

    # bins should tile total chrom length
    assert bins[-1] - 2e4 <= chromosome_length(chrom)

    # convert from bin start to bin end for each bin
    bins = np.concatenate((bins[1:], [chromosome_length(chrom)]))

    # get the segment length of each bin
    segs = bins - np.concatenate(([0], bins[:-1]))

    # mask sites with too little data
    mask = np.isnan(srate)

    # interpolate rate estimates for masked bins from neighboring bins
    srate[mask] = np.interp(bins[mask], bins[~mask], srate[~mask])

    return segs.astype('u8'), srate


def prepare_umap_3(chrom, sub_counts):
    """mini-function to prepare the umap segments and values"""
    bins, bases, mismatch = sub_counts.astype('f8').T
    # highest bin start should be less than chrom length
    assert bins[-1] < chromosome_length(chrom)

    # convert from bin start to bin end for each bin
    bins = np.concatenate((bins[1:], [chromosome_length(chrom)]))
    # get the segment length of each bin
    segs = bins - np.concatenate(([0], bins[:-1]))

    # mask sites with too little data
    mask = ~np.isnan(bases)

    # estimate rates in the remaining windows
    rates = np.zeros(shape=len(bins))
    rates[mask] = mismatch[mask] / bases[mask]

    # interpolate rate estimates for masked bins from neighboring bins
    rates[~mask] = np.interp(bins[~mask], bins[mask], rates[mask])

    return segs.astype('u8'), rates


def prepare_umap_4(chrom, sub_counts):
    """mini-function to prepare the umap segments and values"""
    start, end, bases, mismatch = sub_counts.astype('f8').T
    # segments should cover the entire chromosome
    assert start[0] == 0
    assert end[-1] == chromosome_length(chrom)

    # get the segment length of each bin
    segs = end - start

    # there should be no nan segments
    assert not np.any(np.isnan(bases))

    # estimate rates across windows
    rates = mismatch / bases

    return segs.astype('u8'), rates


def combine_segments(chrom, segments, values):
    """
    Align a list of segment arrays that sum to the same length and their
    corresponding value arrays, e.g., s1, s2, s3 = 3 different segment
    arrays, v1, v2, v3 = the associated value arrays for s1, s2, s3.

    Collapse segment arrays into a single segment array reflecting the min
    value switchpoints among the arrays using the following algorithm:

        1) find the segment array with the shortest segment at index 0, smin,
        and record it.
        2) move to the next index in in the segment array that yielded smin,
        keep others are 0
        3) record all values at the current indices of each segment array
        4) subtract smin from the remaining segments at index 0
        5) return to step 1 until each segment of each map has been processed

    segments
    :param chrom: chromosome used as a key to get expected segment length
    :param segments: a list of arrays containing segments summing to chlen
    :param values: a list of arrays of values corresponding to segmentsn
    :return compressed_map: a map delineated by minimal segmenent breaks
    and including a row of all input map values for each new segment
    """
    # create an array of indices for each data segment
    idx = np.zeros(shape=len(segments), dtype='u8')

    # record each time an index is updated
    steps = 0

    # find min uniform segments of BS/CS and nu values
    newsegs, newvals = [], []
    while any(i < len(s) - 1 for (i, s) in izip(idx, segments)):
        # save min length segment in current row and with current values
        smin = min(s[i] for (i, s) in izip(idx, segments))
        newsegs.append(smin)

        current_values = []
        for (n, i) in enumerate(idx):
            # save current values from each map
            current_values.append(values[n][i])
            # subtract the minimum length from the row of segments
            segments[n][i] -= smin
            # update the index for segments that are reduced to 0
            if segments[n][i] == 0:
                idx[n] += 1

        # save list of current values for the min segment
        newvals.append(current_values)

        # assert at least one step has been taken; update step counter
        assert idx.sum() > steps
        steps = idx.sum()

    # take the final segment lengths and values
    smin = min(s[i] for (i, s) in izip(idx, segments))
    assert all(s[i] == smin for (i, s) in izip(idx, segments))

    # append the fine segment lengths and values
    newsegs.append(smin)
    newsegs = np.array(newsegs)
    assert newsegs.sum() == chromosome_length(chrom)

    newvals.append([v[i] for (i, v) in izip(idx, values)])
    newvals = np.array(newvals)

    return newsegs, newvals


def join_pop_snps(snp_1, snp_2):
    # combine all positions into single array
    all_pos = np.unique(np.concatenate((snp_1[:,0], snp_2[:,0])))

    # get the sample size for filling in gaps
    sample_1 = np.sum(snp_1[0, 1:])
    sample_2 = np.sum(snp_2[0, 1:])
    print sample_1, sample_2

    # get indices of existing SNPs for each pop
    si_1 = np.in1d(all_pos, snp_1[:,0])
    si_2 = np.in1d(all_pos, snp_2[:,0])

    # refill pop 1 arrays with existing poly and mono sites
    s_1 = np.zeros(shape=(all_pos.size, 3))
    # fill in all positions
    s_1[:,0] = all_pos
    # fill in existing SNPs
    s_1[:,1:][si_1] = snp_1[:,1:]
    # fill in monomorphic sites with sample size
    s_1[:,1:][~si_1] = [sample_1, 0]

    # refill pop 2 arrays with existing poly and mono sites
    s_2 = np.zeros(shape=(all_pos.size, 3))
    # fill in all positions
    s_2[:,0] = all_pos
    # fill in existing SNPs
    s_2[:,1:][si_2] = snp_2[:,1:]
    # fill in monomorphic sites with sample size
    s_2[:,1:][~si_2] = [sample_2, 0]

    return s_1, s_2


def compress_data(neutmask, snpcnt, segments):
    """
    compresses polymorphism data array by grouping sites into discrete
    prediction windows
    """
    # create neutral SNPs array from neutmask
    snps = neutral_snps(neutmask, snpcnt)

    # check that every position in SNPs has the same sample size
    sample = np.sum(snps[:, 1:], axis=1)
    assert sample.min() == sample.max()

    # use this for calculations below
    sample = sample[0]
    combos = 0.5 * (sample**2 - sample)

    # use 2 to indicate SNPs in neutmask
    nmsk = neutmask.astype(int)
    nmsk[snps[:,0]] = 2
    assert np.sum(neutmask) == (np.sum(nmsk) - len(snps))

    # create start/end coords for segments
    end = np.cumsum(segments).astype('u4')
    start = np.concatenate(([0], end[:-1])).astype('u4')

    # array for data summaries in each segment
    summaries = np.zeros(shape=(len(end), 5), dtype='f8')

    # go through each segment and get the neutral and poly data for that seg
    snp_idx = 0
    seen_snps = 0
    for i in xrange(len(end)):
        st, en = start[i], end[i]
        # get the neutmask block for current segment
        nblk = nmsk[st:en]
        # get the match and mismatch sites (match just = neut sites, mismatch=0)
        match, mismatch = np.sum(nblk!=0), 0.0
        total = match + mismatch
        # get the number of SNPs in the block
        num_snps = np.sum(nblk==2)
        snp_jdx = snp_idx + num_snps
        # polymorph is the number of SNPs in the segment
        polymorph = snp_jdx - snp_idx
        seen_snps += polymorph
        if polymorph:
            alt = snps[snp_idx:snp_jdx, -1].T
            hom, het = map(np.sum, count_all_pairs(sample, alt))
        else:
            hom = het = 0.0
        # increment hom pairs for all viable sites less poly sites
        if total:
            hom += (combos * (total - polymorph))

        # record sum of homs, hets, viable sites and outgroup mismatches
        pairs = hom + het
        assert (pairs / combos == total) and (pairs % combos == 0)
        summaries[i] = hom, het, total, mismatch, polymorph
        # update the snp index
        snp_idx = snp_jdx

    # check that all snps were added
    assert seen_snps == len(snps)
    # check that neutral site counts in segments sum to all neutral sites
    assert np.sum(summaries[:,2]) == np.sum(neutmask)

    # # convert segment lengths to a set of bins
    # segments = np.concatenate(([0], np.cumsum(segments)))
    # # use histogram on SNP positions to get SNP site indices per seg idx
    # numbins = np.maximum(0, segments-1)
    # # numbins = segments
    # sidx, nidx = np.histogram(a=snps[:, 0], bins=numbins)
    # sidx = np.cumsum(np.concatenate(([0], sidx)))
    # n = len(sidx) - 1  # number of segments
    #
    # # array for data summaries in each segment
    # summaries = np.zeros(shape=(n, 5), dtype='f8')
    #
    # for i in xrange(n):
    #     # indices of snps in current segment of the chrom
    #     si, sj = sidx[i:i+2]
    #     # polymorph = the number of snps in the segment
    #     polymorph = sj - si
    #     # the snp indices pull data from "snps" for the segment
    #     if polymorph:
    #         # compute hom/het pairings possible for the set of snps (if any)
    #         alt = snps[si:sj, -1].T
    #         hom, het = map(np.sum, count_all_pairs(sample, alt))
    #     else:
    #         # if the segment does not contain any snps it is monomorphic
    #         hom = het = 0.0
    #
    #     # neutral mask coordinates of the current segment
    #     ni, nj = nidx[i:i+2].astype(int)
    #     # count match (1+3) and mismatch (2+6) in neutmask block
    #     nblk = neutmask[ni:nj]
    #     match = np.sum((nblk == 1) | (nblk == 3))
    #     mismatch = np.sum((nblk == 2) | (nblk == 6))
    #     # count neutral outgroup matches=1 & mismatches=2 on current segment
    #     # match, mismatch = np.bincount(neutmask[ni:nj], minlength=3)[1:]
    #     total = match + mismatch
    #
    #     # increment hom pairs for all viable sites less poly sites
    #     if total:
    #         hom += (combos * (total - polymorph))
    #
    #     # record sum of homs, hets, viable sites and outgroup mismatches
    #     pairs = hom + het
    #     assert (pairs / combos == total) and (pairs % combos == 0)
    #     summaries[i] = hom, het, total, mismatch, polymorph
    #
    # # print np.sum(summaries[:,2]), np.sum(neutmask)
    # assert np.sum(summaries[:,2]) == np.sum(neutmask)
    # # assert False

    return summaries.T, sample

