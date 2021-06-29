__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv
from itertools import izip
from likelihood.cllh_functions import predicted_pi
from classes.runstruct import ChromStruct, root_dir, cst_from_fldr
from precalc.lh_inputs import load_saved, adjust_arrays
from data_processing.functions import rsquared_function, rsquared_residuals


def standardwindows():
    # use less dense range of scales
    fmcv = '{}/result/rsq_maps/old_data_rsq/phylo-2/' \
           'mcvicker_map_BS1_CS0_161129211053_rsqmap.txt'.format(root_dir)
    return np.loadtxt(fmcv)[:, 0]


def calc_rsquared(scale, cum_pos, ms, nt, dv, pred, residuals=False):
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    obs, prd, num = [], [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate average pi for the window
            pi = nt[i:j, 1].sum() / nt[i:j].sum()
            # calculate average predicted pi for window (weighted by sites)
            pr = np.average(pred[i:j], weights=dv[i:j])
            # calculate site count for window
            n = dv[i:j].sum()
            # record stats to lists
            obs.append(pi)
            prd.append(pr)
            num.append(n)

    # convert all results to arrays
    obs = np.array(obs)
    prd = np.array(prd)
    num = np.array(num)
    msk = (num > 0.2*np.mean(num))  # TODO: what to do about this filter??

    if residuals:
        # calculate residuals
        resd = rsquared_residuals(obs[msk], prd[msk])
        return resd
    else:
        # calculate R squared
        rsq = rsquared_function(obs[msk], prd[msk])
        return rsq


def get_arrays(scale, cum_pos, ms, nt, dv, pred):
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    obs, prd, num = [], [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate average pi for the window
            pi = nt[i:j, 1].sum() / nt[i:j].sum()
            # calculate average predicted pi for window (weighted by sites)
            pr = np.average(pred[i:j], weights=dv[i:j])
            # calculate site count for window
            n = dv[i:j].sum()
            # record stats to lists
            obs.append(pi)
            prd.append(pr)
            num.append(n)

    # convert all results to arrays
    obs = np.array(obs)
    prd = np.array(prd)
    num = np.array(num)
    msk = (num > 0.2*np.mean(num))  # TODO: what to do about this filter??

    return np.column_stack((prd[msk], obs[msk]))


def get_yri_params(fldr):
    """get just the predictions for YRI map"""
    cst = cst_from_fldr(fldr)
    return cst.params, cst.stat.meanpi


def main():
    if len(argv) != 2:
        print 'usage: calc_rsq <folder_name>'
        exit(1)

    # set folder path and get ChromStruct from folder
    fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
    cst = cst_from_fldr(argv[1])

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

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

    if 'cadd' in argv[1]:
        # yri_fldr = 'cadd94_gmask_mnb_378'
        yri_fldr = 'cadd94_gmask_v1.6_without_bstat'
    else:
        # assert 'fish_cons' in argv[1]
        yri_fldr = 'fish_cons94_new'

    # get bmap params from YRI, use theta from specific population
    yri_prm, yri_meanpi = get_yri_params(yri_fldr)

    # get prediction from YRI map
    yri_pred = predicted_pi(yri_prm, cst, nu, bs, cs)
    # rescale predictions to fit the mean pi of current population
    yri_pred *= (cst.stat.meanpi / yri_meanpi)
    # get prediction from current population params
    pred = predicted_pi(cst.params, cst, nu, bs, cs)

    # # use init file path as template for rsq file
    # rsq_file = fdir +'/rsq.YRI.map.log'
    # # rsq_file = fdir +'/rsq.CEU.map.log'
    # residuals_file = fdir + '/YRI.residuals.log'

    # # calculate r squared
    # with open(rsq_file, 'w') as f:
    #     for sc in standardwindows():
    #         rsq = calc_rsquared(sc, cum_pos, msk, nt, dv, pred)
    #         f.write('{:.3e}\t{}\n'.format(sc, rsq))

    # with open(residuals_file, 'w') as f:
    #     for sc in standardwindows():
    #         rsq = calc_rsquared(sc, cum_pos, msk, nt, dv, pred, residuals=True)
    #         f.write('{:.3e}\t{}\n'.format(sc, rsq))

    # array files
    if not os.path.isdir(fdir + 'arrays'):
        os.mkdir(fdir + 'arrays')
    for sc in standardwindows():
        # first array file is for YRI map result
        array_file_1 = fdir + 'arrays/yri_pred_obs_{:.3e}.npz'.format(sc)
        prd_obs_1 = get_arrays(sc, cum_pos, msk, nt, dv, yri_pred)
        np.savez_compressed(array_file_1, prdobs=prd_obs_1)
        # second array file is for map from current population params
        array_file_2 = fdir + 'arrays/pred_obs_{:.3e}.npz'.format(sc)
        prd_obs_2 = get_arrays(sc, cum_pos, msk, nt, dv, pred)
        np.savez_compressed(array_file_2, prdobs=prd_obs_2)


if __name__ == '__main__':
    main()
