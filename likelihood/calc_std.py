__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv
from itertools import izip
from likelihood.cllh_functions import predicted_pi
from classes.runstruct import ChromStruct, root_dir
from precalc.lh_inputs import load_saved, adjust_arrays
from data_processing.functions import rsquared_function, swap_root


def standardwindows():
    # use less dense range of scales
    fmcv = '{}/result/rsq_maps/old_data_rsq/phylo-2/' \
           'mcvicker_map_BS1_CS0_161129211053_rsqmap.txt'.format(root_dir)
    return np.loadtxt(fmcv)[:, 0]


def calc_std(cum_pos, ms, nt, dv):
    scale = 1e6
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    obs, num = [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate average pi for the window
            pi = nt[i:j, 1].sum() / nt[i:j].sum()
            # calculate site count for window
            n = dv[i:j].sum()
            # record stats to lists
            obs.append(pi)
            num.append(n)

    # convert all results to arrays
    obs = np.array(obs)
    num = np.array(num)
    msk = (num > 0.2*np.mean(num))  # TODO: what to do about this filter??

    # calculate standard deviation
    std = np.std(obs[msk])

    return std


def calc_bvar(cum_pos, ms, bsx, dv):
    """get the variance in B values at 1Mb scale"""
    scale = 1e6
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    bval, num = [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate mean B for window
            bmean = np.average(bsx[i:j], weights=dv[i:j])
            # calculate site count for window
            n = dv[i:j].sum()
            # record stats to lists
            bval.append(bmean)
            num.append(n)

    # convert all results to arrays
    bval = np.array(bval)
    num = np.array(num)
    msk = (num > 0.2*np.mean(num))  # TODO: what to do about this filter??

    # calculate standard deviation
    bvar = np.var(bval[msk])

    return bvar


def calc_predvar(cum_pos, ms, pred, dv):
    """get the variance in B values at 1Mb scale"""
    scale = 1e6
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    predval, num = [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate mean B for window
            predmean = np.average(pred[i:j], weights=dv[i:j])
            # calculate site count for window
            n = dv[i:j].sum()
            # record stats to lists
            predval.append(predmean)
            num.append(n)

    # convert all results to arrays
    predval = np.array(predval)
    num = np.array(num)
    msk = (num > 0.2*np.mean(num))  # TODO: what to do about this filter??

    # calculate standard deviation
    predvar = np.var(predval[msk])

    return predvar


def main():
    if len(argv) != 2:
        print 'usage: calc_std <folder_name>'
        exit(1)

    fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
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
    # print msk.shape, mnu1.shape
    msk &= mnu1[:,0]

    # get bmap to get percentiles
    params = cst.params
    uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
    bwt = uvec / cst.fixed.u_fix
    bsx = np.exp(np.dot(bs, bwt))

    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(cst.params, cst, nu, bs, cs)

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

    # use init file path as template for rsq file
    std_file = '/'.join(f_init.split('/')[:-1]) +'/std.log'
    var_file = '/'.join(f_init.split('/')[:-1]) +'/bvar.log'
    predvar_file = '/'.join(f_init.split('/')[:-1]) + '/predvar.log'

    # calculate std dev
    with open(std_file, 'w') as f:
        std = calc_std(cum_pos, msk, nt, dv)
        f.write('{:.3e}\t{}\n'.format(1e6, std/cst.stat.meanpi))

    # calculate b variance
    with open(predvar_file, 'w') as f:
        # bvar = calc_bvar(cum_pos, msk, bsx, dv)
        predvar = calc_predvar(cum_pos, msk, pred, dv)
        f.write('{:.3e}\t{}\n'.format(1e6, predvar / pred.mean()**2))


if __name__ == '__main__':
    main()
