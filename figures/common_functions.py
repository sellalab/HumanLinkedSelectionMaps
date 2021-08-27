__author__ = 'davidmurphy'

import numpy as np
from skmisc import loess
from data_processing.functions import rsquared_function
from classes.runstruct import root_dir, chromosome_length, cst_from_fldr


final_dir = root_dir + '/result/final_files'


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def get_bbins(tkn=None, nbins=100):
    if nbins != 100:
        sdir = root_dir + '/compress/sorted_bins/sortbin{}'.format(nbins)
    else:
        sdir = root_dir + '/compress/sorted_bins'

    bins = []
    for i in xrange(nbins):
        if tkn:
            f_name = sdir + '/sorted_bin{}.{}.npz'.format(i, tkn)
        else:
            f_name = sdir + '/sorted_bin{}.npz'.format(i)
        bins.append(np.load(f_name)['sg'])

    return bins


def format_panels(ax):
    # axis format variables for white background
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25
    ax.set_facecolor('white')
    ax.grid(**gridargs)
    ax.axis('on')
    ax.patch.set_edgecolor('black')
    ax.patch.set_linewidth(axis_lw)


def rsq_three_scales(fldr):
    r_file = final_dir + '/{}/rsq.log'.format(fldr)
    w, r = np.loadtxt(r_file).T
    ikeep = np.array([0, 7, 13])

    return r[ikeep]


def get_sortplot_rsq(fldr, num):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    sort_file = fdir + 'basic_sort_n{}.txt'.format(num)
    div, pi, pred = np.loadtxt(sort_file).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    # normalize by pi0
    pi /= pi0
    pred /= pi0

    return rsquared_function(pi, pred)


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


def hist_error(eps, nbins, rlist):
    """function used to create histogram plots of map error"""
    tol = 2.0 * eps  # tolerance range = 2 x epsilon
    bnd = 2.0 * tol  # plot bounds = 4 x epsilon
    inc = 2.0 * bnd / nbins  # bin increments according to number of bins
    bns = np.arange(-bnd, bnd+inc, inc)

    hlist = []  # list for histogram counts of error values
    olist = []  # list for outliers
    out = tol * 2.0 ** np.arange(5)
    for r in rlist:
        # make a histogram of errors within +/-4 x epsilon
        histo = np.histogram(r, bns)
        hlist.append(histo)
        # count outliers above log range of thresholds
        outlr = np.array([100.0 * np.sum(abs(r) > o) / r.size for o in out])
        outlr[:-1] -= outlr[1:]
        olist.append(outlr)

    return hlist, olist


def get_pop_array(win, pop, anno, arr_type='YRI'):
    """get obs/pred arrays from self-map or YRI-map for a single pop"""
    # set file name depending on array type being used
    if arr_type == 'YRI':
        f_str = '/arrays/yri_pred_obs_{:.3e}.npz'.format(win)
    else:
        f_str = '/arrays/pred_obs_{:.3e}.npz'.format(win)

    if pop == 'YRI':
        # if ('94' in anno) and ('without_bstat' not in anno):
        #     d_str = '/{}_mnb_378'.format(anno)
        # else:
        #     d_str = '/{}'.format(anno)
        d_str = '/{}'.format(anno)
        fa = final_dir + d_str + f_str
    else:
        d_str = '/{}_{}'.format(pop, anno)
        fa = final_dir + d_str + f_str

    prd, obs = np.load(fa)['prdobs'].T

    return obs, prd


def get_pop_array_2(win, fldr, arr_type='YRI'):
    """get obs/pred arrays from self-map or YRI-map for a single pop"""
    # set file name depending on array type being used
    if arr_type == 'YRI':
        f_str = '/arrays/yri_pred_obs_{:.3e}.npz'.format(win)
    else:
        f_str = '/arrays/pred_obs_{:.3e}.npz'.format(win)

    fa = final_dir + '/' + fldr + f_str

    prd, obs = np.load(fa)['prdobs'].T

    return obs, prd


def get_centromere_dict():
    centro_dict = {}
    gfile = root_dir + '/data/coords/gaps/hg19_genome_gaps.txt'
    with open(gfile, 'r') as f:
       f.next()
       for line in f:
           l = line.split()
           ch, start, end, name = l[1], int(l[2]), int(l[3]), l[-2]
           if name == 'centromere':
                centro_dict[ch] = (start, end)
                relstart = float(start) / chromosome_length(ch)
                relend = float(end) / chromosome_length(ch)
                print '{} {} {}'.format(ch, relstart, relend)

    return centro_dict


def get_telomere_dict():
    telo_dict = {}
    for c in xrange(1, 23):
        # set initial max to 0, initial min to 3e8
        ch = 'chr{}'.format(c)
        telo_dict[ch] = (1e4, chromosome_length(ch)-1e4)

    return telo_dict

