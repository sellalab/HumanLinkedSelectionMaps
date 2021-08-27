from classes.runstruct import *
from obsolete.calculate_likelihood import calc
from obsolete.preprocess import pre_process_arrays

__author__ = 'davidmurphy'


"""
A collection of functions for analyzing genomic maps of diversity (predicted and observed)
"""


def binned_map(vectors, binsize, stepsize=None, normalize=False, weights=None):
    """
    This function uses the params stored in a RunStruct to construct a map for the chromosome specified and then
    calculates the mean observed and predicted diversity levels in that map in binsize windows. If a step is
    specified, the means are calculated in sliding windows of size binsize with 'step' degree of overlap.
    :param vectors: tuple or list of vectors returned by inferred_map function -- pos, hom, het, u, unweighted_pred
    :param binsize: the number of bp that should be used to average over
    :param stepsize: the overlap for sliding windows (default=none)
    :param normalize: if True, normalize the vector of observed and predicted diversity levels to 1 (default=False)
    :param weights: a vector of weights for each data point (default is none -- equal weighting)
    :return pos, obs, pred: the upper edge of each binsize windows used to tile the chromosome, observed and
    predicted mean diversity levels in each window
    """
    # upper edge positions, hom/het pair counts, u-rate and unweighted prediction for each segment with neutral sites
    pos, hom, het, u, unweighted_pred = vectors

    # apply weights to the hom/het counts
    if weights is not None:
        het *= weights
        hom *= weights

    # scale prediction by pairs in each segment - predictions in site-rich segments given greater weight
    pairs = het + hom  # total possible pairs per segment
    weighted_pred = pairs * unweighted_pred
    # get the sum of pairs, hets and sample weighted prediction in each bin
    bin_pos, bin_pairs = binned_function(pos, pairs, binsize, np.sum, stepsize)
    _, bin_hets = binned_function(pos, het, binsize, np.sum, stepsize)
    _, bin_pred = binned_function(pos, weighted_pred, binsize, np.sum, stepsize)
    low_idx = (bin_pairs < 0.05 * np.nanmean(bin_pairs))
    # NOTE: just turn aberrant pairs to nans, smoothen does not seem to work in all situations...
    bin_pairs[low_idx] = np.nan
    # turning 0's to nans improves the rsquare for some reason (and the look of the plot)
    bin_hets[bin_hets == 0] = np.nan
    # correct for aberrant bins where the number of pairs is far below the mean by merging them with adjacent bins
    # if low_idx.any():
    #     smooth_edges(bin_pairs, bin_hets, bin_pred, low_idx)
    # mean pi for a bin is same as for a site: total het pairs / total pairs
    mean_pi = bin_hets / bin_pairs
    # mean prediction is a weighted average of all segments in the bin
    mean_pred = bin_pred / bin_pairs
    # optional normalize
    if normalize:
        # normalize mean pi and prediction values to 1 by dividing them with their respective means across bins
        return bin_pos, mean_pi / np.nanmean(mean_pi), mean_pred / np.nanmean(mean_pred)
    else:
        return bin_pos, mean_pi, mean_pred


def smooth_edges(pairs, hets, weighted_pred, low_idx):
    """
    Function that merges bins with very little data into adjacent data to avoid huge swings in calculated rsquared
    values that such bins can lead to
    :param pairs: total possible pairs in the bin
    :param hets: total het pairs in bin
    :param weighted_pred: prediction for the bin weights by pairs
    :param low_idx: the positions of bins with very little data
    """
    low_idx = np.where(low_idx)[0]
    # put half of the data in the previous bin and the other half in the following bin, turn low data bins to nan
    for i in low_idx:
        # indices must not extend beyond the array length
        ilo = max(0, i - 1)
        ihi = min(len(pairs) - 1, i + 1)
        for data in [pairs, hets, weighted_pred]:
            data[ilo] += (0.5 * data[i])
            data[ihi] += (0.5 * data[i])
            data[i] = np.nan
    return None


def binned_function(positions, values, binsize, function=None, step=None):
    """
    A function to run some function on 'values' in bins of size 'binsize' using 'positions' to determine which bin a
    value should be placed in. The default function is numpy.mean. Optional 'step' can be used for overlapping bins,
    where the default step=binsize with no overlap between bins.
    :param positions: a vector of positions or ranks in the same units as 'binsize'
    :param values: a set of numerical values collected in bins (should correspond 1:1 with positions) on which
    'function' is applied.
    :param binsize: the size of the bins to use
    :param function: a function to use on the binned data, which will be a 1D numpy.ndarray. default is numpy.mean
    :param step: an optional parameter to step from bin_n to bin_n+1 in some increment that is < binsize (for rolling
    averages for example).
    :return edges, f_values: the upper bounds corresponding to each bin and the function value for that bin
    """

    # start the bottom edge of the bins at min(positions) and continue until the upper edge of the final bin is >=
    # max(positions): numpy.histogram will include points in the final bin that are equal to the LAST edge in the set
    # of edges, although all other bins include values where: edge_lower <= value < edge_upper
    bmin, bmax = min(positions), max(positions)
    edges = np.arange(start=bmin, stop=bmax + binsize, step=binsize)

    counts = []
    bounds = []
    if step:
        # for this formulation binsize must be evenly divisible by step. step is adjusted up or down to make sure
        num_steps = int(0.5 + binsize / step)
        step = 1.0 * binsize / num_steps
        shift = 0
        for i in xrange(num_steps):
            n, b = np.histogram(a=positions, bins=edges + shift)
            counts.append(n.astype('f8'))
            bounds.append(b[1:])
            shift += step
    else:
        # count the number of 'positions' in each bin, use to index 'values' by bin
        n, b = np.histogram(a=positions, bins=edges)
        counts.append(n)
        bounds.append(b[1:])

    start = 0
    f_values = []
    for n in counts:
        for i in xrange(len(n)):
            end = start + n[i]
            # filter for 'nan'
            if np.isfinite(values[start:end]).any():
                mask = np.isfinite(values[start:end])
                # use input function if specified or default to mean
                if function:
                    f_x = function(values[start:end][mask])
                else:
                    f_x = np.mean(values[start:end][mask])
                f_values.append(f_x)
            # if all 'nan', record 'nan' for the average in that bin
            else:
                f_values.append(np.nan)
            # start at the next iteration from the end of the current iteration
            start = end

        # set after the first step in the previous grid for starting point of the next grid
        start = n[0]

    bounds = np.concatenate(bounds)
    f_values = np.array(f_values, dtype='f8')
    idx = bounds.argsort()
    # if there was just 1 set of edges, don't take the first value
    # if not step:
    #     edges = edges[1:]

    return bounds[idx], f_values[idx]


def infmap_chrom(rst, chrom=None, uconst=False, div=False, con=False, gc=False, bt=False):
    """
    A function to create a final map of **predicted diversity** by applying the params from the inference to the array
    of B-maps to produce a single B-map scaled from 0-1 in each bin. For each B value in each bin j we calculated
    predicted diversity by multiplying Bj * pi0 * mu_local.
    :param rst: RunStruct from the completed run, which contains the inferred params, file directories, etc.
    :param chrom: the chromosome for which to generate the map.
    :param uconst: a flag to set if a constant mutation rate correction should be used for predictions only
    :param div: a flag to return divergence with pre_process_arrays
    :param con: a flag to return conservation scores with pre_process_arrays
    :param gc: a flag to return gc-content with pre_process_arrays
    :param bt: a flag to return the bootstrap sample weights for each data point
    :return inferred_map: an array in McVicker-type coordinates with predicted reduction factor (0 to 1 scale) per seg
    """
    assert isinstance(rst, RunStruct)
    # turn on result mode in case it is not already on for the RunStruct
    rst.vars.complete = True
    # use chrom to specify the file indices to use
    chidx = rst.chroms.index(chrom)
    # configure the options for pre-processing
    options = dict(return_pos=1, return_rst=0, return_div=div, return_con=con, return_gc=gc, return_bootstrap=bt)
    # pre-process arrays for the current chromosome, retrieve positions
    arrs = pre_process_arrays(rst, chrom_range=xrange(chidx, chidx + 1), **options)
    nt, u, bs, cs, pos = arrs[:5]
    hom, het = nt.T

    # if the uconst flag is applied, create a vector of meandiv to replace the vector u for a constant mutation rate
    if uconst:
        umean = np.ones(shape=len(u)) * rst.stat.meandiv
        prediction = calc(u=umean, hom=hom, het=het, rst=rst, bsx=bs, csx=cs)
    else:
        prediction = calc(u=u, hom=hom, het=het, rst=rst, bsx=bs, csx=cs)

    # return positions, poly data and predictions together, as well as any optional arrays
    return tuple([pos, hom, het, u, prediction] + list(arrs[5:]))


def infmap_full(rst, uconst=False, div=False, con=False, gc=False, bt=False):
    """
    Return the full map and any optional data indicated
    :param rst: RunStruct from the completed run, which contains the inferred params, file directories, etc.
    :param uconst: a flag to set if a constant mutation rate correction should be used for predictions only
    :param div: a flag to return divergence with pre_process_arrays
    :param con: a flag to return conservation scores with pre_process_arrays
    :param gc: a flag to return gc-content with pre_process_arrays
    :param bt: a flag to return the bootstrap sample weights for each data point
    :return inferred_map: an array in McVicker-type coordinates with predicted reduction factor (0 to 1 scale) per seg
    """
    assert isinstance(rst, RunStruct)
    # turn on result mode in case it is not already on for the RunStruct
    rst.vars.complete = True
    # configure the options for pre-processing
    options = dict(return_bootstrap=bt, return_div=div, return_con=con, return_gc=gc)
    arrs = pre_process_arrays(rst, **options)
    nt, u, bs, cs = arrs[:4]  # extract basic arrays
    hom, het = nt.T
    # get the bootstrap out if it is flagged
    if bt:
        s = arrs[4]
    else:
        s = None
    # if the uconst flag is applied, create a vector of meandiv to replace the vector u for a constant mutation rate
    if uconst:
        umean = np.ones(shape=len(u)) * rst.stat.meandiv
        prediction = calc(u=umean, hom=hom, het=het, rst=rst, bsx=bs, csx=cs, s=s)
    else:
        prediction = calc(u=u, hom=hom, het=het, rst=rst, bsx=bs, csx=cs, s=s)

    # full data set without positions
    return tuple([hom, het, u, prediction] + list(arrs[4:]))


def infmap_stats(rst, bootstrap=False):
    """
    Calculate statistics on the completed map
    :param rst: RunStruct
    :param bootstrap: use bootstrap flag
    """
    hom, het, u, pred, boot = infmap_full(rst=rst, bt=bootstrap)
    bootmask = boot > 0
    rst.stat.meanboot = np.mean(boot[bootmask])
    rst.stat.varboot = np.var(boot[bootmask])
    rst.stat.meanpred = np.sum(boot * pred) / np.sum(boot)
    rst.stat.varpred = np.sum(boot * ((pred - rst.stat.meanpred) ** 2)) / np.sum(boot)
    rst.calc_stats()
    rst.save()

    return None


def bin_bmap_scores(rst):
    assert isinstance(rst, RunStruct)
    # np.log(np.maximum(b, 1) / 100.0))
    barr = np.concatenate([np.log(np.maximum(np.load(f), 1) / 100.0) for f in rst.bs_files])
    sarr = np.concatenate([np.load(f) for f in rst.seg_files])
    bwt = np.power(10, rst.params[rst.fixed.bi:rst.fixed.bj]) / rst.fixed.u_fix
    bs = np.exp(np.dot(barr, bwt))
    bs *= 100.0
    cnts, bins = np.histogram(a=bs, bins=np.arange(0, 100.0, 1.0), weights=sarr)
    return np.column_stack([bins[1:], cnts])
