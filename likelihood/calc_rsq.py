__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv
from likelihood.cllh_functions import predicted_pi
from classes.runstruct import ChromStruct, root_dir
from precalc.lh_inputs import load_saved, adjust_arrays
from data_processing.functions import rsquared_function, swap_root


def standardwindows():
    # use less dense range of scales
    fmcv = '{}/result/rsq_maps/old_data_rsq/phylo-2/' \
           'mcvicker_map_BS1_CS0_161129211053_rsqmap.txt'.format(root_dir)
    return np.loadtxt(fmcv)[:, 0]


def calc_rsquared(scale, cum_pos, ms, nt, dv, pred):
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    obs, prd, num = [], [], []
    for (i, j) in zip(idx, jdx):
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

    # calculate R squared
    rsq = rsquared_function(obs[msk], prd[msk])

    return rsq


def main():
    if len(argv) != 2:
        print('usage: calc_rsq <folder_name>')
        exit(1)
    # folder_name = argv[1]
    # fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
    # # list of result files
    # flst = [fdir + f for f in os.listdir(fdir) if f.endswith('.txt')]
    # # list of chromstructs for each result
    # rlst = [ChromStruct('chr1', init=f) for f in flst]
    #
    # # create tuples of file index, final LH result
    # ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
    # # get indices of the top 3 best LH runs (after sorting on LH)
    # top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]
    # # take the average of the params of the top 3 runs
    # pavg = np.average([rlst[i].params for i in top3], axis=0)
    #
    # # use any file for init
    # init = flst[0]
    # swap_root(init)
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
    nu, nt, dv, pl = [a[mnu2] for a in [nu, nt, dv, pl]]
    # print msk.shape, mnu1.shape
    msk &= mnu1[:,0]

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(cst.params, cst, nu, bs, cs)
    # print 'max pred: {}'.format(pred.max())
    # exit(1)

    # apply error rate c threshold if in use
    if cst.fixed.cth is not None:
        c = cst.fixed.cth
        # get the standard prediction
        pred = predicted_pi(cst.params, cst, nu, bs, None)
        # get divergence adjusted mean pi
        pii_avg = nu * (cst.stat.meanpi / cst.stat.meandiv)
        # change prediction to error weighted prediction
        pred = ((1.0 - c) * pred) + (c * pii_avg)

    # use init file path as template for rsq file
    rsq_file = '/'.join(f_init.split('/')[:-1]) +'/rsq.log'
    # calculate r squared
    with open(rsq_file, 'w') as f:
        for sc in standardwindows():
            rsq = calc_rsquared(sc, cum_pos, msk, nt, dv, pred)
            f.write('{:.3e}\t{}\n'.format(sc, rsq))


if __name__ == '__main__':
    main()
