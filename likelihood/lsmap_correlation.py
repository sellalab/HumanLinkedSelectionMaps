__author__ = 'davidmurphy'

import os
import numpy as np
from itertools import izip
from scipy.stats import pearsonr
from sys import argv, stdout, stderr
from precalc.lh_inputs import load_saved, adjust_arrays
from classes.runstruct import ChromStruct, root_dir
from likelihood.cllh_functions import predicted_pi
from likelihood.calc_rsq import standardwindows


def binned_pred(scale, cum_pos, ms, dv, pred):
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # calculate mean predicted value in bins
    prd = []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate average predicted pi for window (weighted by # sites)
            pr = np.average(pred[i:j], weights=dv[i:j])
            prd.append(pr)
        else:
            prd.append(0)

    prd = np.array(prd)

    return prd


def get_pred(fldr, scale):
    """get map of predictions for a given results folder"""
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)
    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

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
    msk &= mnu1[:, 0]

    # use constant mutation rate for predictions!
    nu_pred = np.full(len(nu), cst.stat.meandiv)

    # load prediction map
    pred = predicted_pi(cst.params, cst, nu_pred, bs, cs)

    return binned_pred(scale, cum_pos, msk, dv, pred)


def lsmap_autocor(fldr_1, fldr_2):
    """calculate the Pearson correlation for two sets of predictions"""

    for w in standardwindows():
        pred_1 = get_pred(fldr_1, w)
        pred_2 = get_pred(fldr_2, w)
        imax = min(pred_1.size, pred_2.size)
        pear = pearsonr(pred_1[:imax], pred_2[:imax])[0]
        m = '{:.4e} {} sizes {} {}\n'.format(w, pear, pred_1.size, pred_2.size)
        stderr.write(m)
        stdout.flush()


def main():
    if len(argv) != 3:
        print 'usage: lsmap_correlation <fldr_1> <fldr_2>'
        exit(1)
    lsmap_autocor(argv[1], argv[2])


if __name__ == '__main__':
    main()

