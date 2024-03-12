__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from precalc.lh_inputs import load_saved, adjust_arrays
from classes.runstruct import ChromStruct, root_dir
from likelihood.cllh_functions import predicted_pi


def basic_sort(nt, ns, nu, pred, num):
    """Sort data arrays by predictions and average results into bins"""
    # get sorting indices based on prediction value
    sidx = np.argsort(pred)
    # sort predictions and div/poly data
    nt, ns, nu, pred = [a[sidx].astype('f8') for a in [nt, ns, nu, pred]]
    # calculate mean div
    meandiv = np.average(nu, weights=ns)

    # get start and end indices per partition
    idx = sortbin_edges(sites=ns, numbins=num)

    # gather data summaries
    sorted_array = []
    for (i, j) in idx:
        # calculate each statistic for the present bin
        # norm_div = np.sum(subs[i:j]) / (meandiv * np.sum(sites[i:j]))
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


def main():
    if len(argv) != 2:
        print('usage: calc_rsq <folder_name>')
        exit(1)

    fldr = argv[1]
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # use leave one out predictions for result
    leave_one_out = False
    dropped_chrom = True

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # # ****** LOAD CEU DATA AS WELL ******
    # ceu_init = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/result/init_files/CEU.ape_cons95.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    # ceu_cst = ChromStruct('chr1', init=ceu_init)
    # ceu_nt = load_saved(ceu_cst)[4]
    # nt = ceu_nt
    # # ******

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
    nu, nt, dv, pl = [a[mnu2] for a in [nu, nt, dv, pl]]
    msk &= mnu1[:,0]

    # use constant mutation rate for predictions!
    nu_pred = np.full(len(nu), cst.stat.meandiv)

    # calculate predicted pi from arrays and average of top 3 params
    # cst.fixed.min_bsx = None
    pred = predicted_pi(cst.params, cst, nu_pred, bs, cs)

    # apply error rate c threshold if in use
    if cst.fixed.cth is not None:
        c = cst.fixed.cth
        # get the standard prediction
        pred = predicted_pi(cst.params, cst, nu_pred, bs, None)
        # get divergence adjusted mean pi
        pii_avg = nu_pred * (cst.stat.meanpi / cst.stat.meandiv)
        # change prediction to error weighted prediction
        pred = ((1.0 - c) * pred) + (c * pii_avg)

    if leave_one_out:
        print('USING LEAVE ONE OUT PREDICTIONS MAP FOR RESULTS FOLDER: {}'.format(fldr))
        f_pred = root_dir + '/result/final_files/{fldr}_jackknife_results/{fldr}.jkpred_const.npy'.format(fldr=fldr)
        pred = np.load(f_pred)[mnu2]
        assert len(pred) == len(nu)

    if dropped_chrom:
        print('USING DROPPED CHROMOSOME PREDICTIONS MAP')
        f_pred = root_dir + '/result/final_files/drop_chrom_results/drop_chrom_pred_nuconst.npy'
        pred = np.load(f_pred)[mnu2]
        assert len(pred) == len(nu)

    # # get bmap to get percentiles
    # params = cst.params
    # min_bsx = cst.fixed.min_bsx
    # uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
    # bwt = uvec / cst.fixed.u_fix
    # bsx = np.exp(np.dot(bs, bwt))
    # # apply a cutoff to the composite bmap
    # if min_bsx:
    #     bsx = np.maximum(bsx, min_bsx)
    # bpct = np.percentile(bsx, np.arange(100))

    # create output filename
    # sort_file = '/'.join(f_init.split('/')[:-1]) + '/predsort_use_bthresh.txt'

    # bpct_file = '/'.join(f_init.split('/')[:-1]) + '/b_percentiles.txt'

    # # ****** SAVE TO CEU DIR ******
    # sort_file = '/'.join(f_init.split('/')[:-2]) + '/ceu_ape95/predsort_use_bthresh.txt'
    # bpct_file = '/'.join(f_init.split('/')[:-2]) + '/ceu_ape95/b_percentiles.txt'
    # # ******

    # get sorted array, default n=100
    # n = 100
    for n in [100, 250, 500, 1000, 2000]:
        sarr = basic_sort(nt, dv, nu, pred, n)
        sort_file = fdir + 'basic_sort_n{}.txt'.format(n)
        if leave_one_out:
            sort_file = fdir + 'basic_sort_n{}_leaveOneOut.txt'.format(n)
        np.savetxt(sort_file, sarr)


if __name__ == '__main__':
    main()
