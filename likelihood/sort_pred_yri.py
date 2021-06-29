__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import load_saved, adjust_arrays
from classes.runstruct import ChromStruct, root_dir, cst_from_fldr


def basic_sort(nt, ns, nu, pred, num):
    """Sort data arrays by predictions and average results into bins"""
    # get sorting indices based on prediction value
    sidx = np.argsort(pred)
    # sort predictions and div/poly data
    nt, ns, nu, pred = [a[sidx].astype('f8') for a in nt, ns, nu, pred]
    # calculate mean div
    meandiv = np.average(nu, weights=ns)

    # get start and end indices per partition
    idx = sortbin_edges(sites=ns, numbins=num)

    # gather data summaries
    sorted_array = []
    for (i, j) in idx:
        # calculate each statistic for the present bin
        norm_div = np.average(nu[i:j], weights=ns[i:j]) / meandiv
        scaled_pi = np.sum(nt[i:j, 1]) / (np.sum(nt[i:j]) * norm_div)
        # use weighted mean for pred
        prediction = np.average(pred[i:j], weights=ns[i:j])
        sorted_array.append([norm_div, scaled_pi, prediction])

    return np.array(sorted_array)


def sortbin_edges(sites, numbins):
    """
    get the upper indices of sorted data array that divide
    data into bins of equal neutral site counts
    """

    # get number of sites needed so that numbins x numsites = total sites
    numsites = int(np.sum(sites) / numbins)

    # find indices that partition cumulative sorted site count into numsites
    cumsites = np.cumsum(sites)
    bounds = np.arange(numsites, cumsites[-1], numsites)

    # get the ending index for each partition
    jdx = list(np.searchsorted(a=cumsites, v=bounds))

    # return a list of (start, end) indices for each partition
    return zip([0] + jdx[:-1], jdx)


def get_yri_params(fldr):
    """get just the predictions for YRI map"""
    cst = cst_from_fldr(fldr)
    return cst.params, cst.stat.meanpi


def main():
    if len(argv) != 2:
        print 'usage: calc_rsq <folder_name>'
        exit(1)

    fldr = argv[1]
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    cutoff = 0.25
    mnu1 = (nu <= cutoff)

    # mask and rescale maps
    bs, cs, nu, nt, dv, pl, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    mnu2 = (nu <= cutoff)
    if bs is not None:
        bs = bs[mnu2]
    if cs is not None:
        cs = cs[mnu2]
    nu, nt, dv, pl = [a[mnu2] for a in nu, nt, dv, pl]
    msk &= mnu1[:,0]

    # use constant mutation rate for predictions!
    nu_pred = np.full(len(nu), cst.stat.meandiv)

    # get prediction params from YRI CADD results
    yri_fldr = 'cadd94_gmask_mnb_378'
    # get bmap params from YRI, use theta from specific population
    yri_prm, yri_meanpi = get_yri_params(yri_fldr)
    # get prediction from YRI map
    pred = predicted_pi(yri_prm, cst, nu_pred, bs, cs)
    # rescale predictions to fit the mean pi of current population
    pred *= (cst.stat.meanpi / yri_meanpi)

    for n in [100, 250, 500, 1000, 2000]:
        sarr = basic_sort(nt, dv, nu, pred, n)
        sort_file = fdir + 'YRI_sorted.basic_sort_n{}.txt'.format(n)
        np.savetxt(sort_file, sarr)


if __name__ == '__main__':
    main()
