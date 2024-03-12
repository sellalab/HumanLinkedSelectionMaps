__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from likelihood.cllh_functions import predicted_pi
from classes.runstruct import ChromStruct, root_dir
from precalc.lh_inputs import load_saved, adjust_arrays


def predicted_and_observed(cum_pos, ms, nt, dv, pred, scale, slide=None):
    # create indices for sliding windows
    if slide:
        assert slide < 1
        windows = np.arange(0, cum_pos[-1], scale * slide)
        si = np.searchsorted(cum_pos[ms], windows)
        idx = si[:-2]
        jdx = si[2:]
    # create indices for nonoverlapping windows
    else:
        # create upper, lower indices to sort data into genomic windows
        windows = np.arange(0, cum_pos[-1], scale)
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
        else:
            obs.append(np.nan)
            prd.append(np.nan)
            num.append(0)

    # convert all results to arrays
    obs = np.array(obs)
    prd = np.array(prd)
    num = np.array(num)

    return np.column_stack((prd, obs, num))


def main():
    if len(argv) != 2:
        print('usage: calc_rsq <folder_name>')
        exit(1)
    scale = 1e6
    slide = 0.5
    chrom = 'chr1'
    leave_one_out = False
    dropped_chrom = True

    fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst, chroms=[chrom])

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
    msk &= mnu1[:, 0]

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)
    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(cst.params, cst, nu, bs, cs)
    if leave_one_out:
        print('USING LEAVE ONE OUT PREDICTIONS MAP FOR RESULTS FOLDER: {}'.format(argv[1]))
        f_pred = root_dir + '/result/final_files/{fldr}_jackknife_results/{fldr}.jkpred.npy'.format(fldr=argv[1])
        pred = np.load(f_pred)[:len(mnu2)][mnu2]
        # print('pred={} nu={}'.format(len(pred), len(nu)))
        assert len(pred) == len(nu)

    if dropped_chrom:
        print('USING DROPPED CHROMOSOME PREDICTIONS MAP')
        f_pred = root_dir + '/result/final_files/drop_chrom_results/drop_chrom_pred.npy'
        pred = np.load(f_pred)[:len(mnu2)][mnu2]
        assert len(pred) == len(nu)

    # apply error rate c threshold if in use
    if cst.fixed.cth is not None:
        c = cst.fixed.cth
        # get the standard prediction
        pred = predicted_pi(cst.params, cst, nu, bs, None)
        # get divergence adjusted mean pi
        pii_avg = nu * (cst.stat.meanpi / cst.stat.meandiv)
        # change prediction to error weighted prediction
        pred = ((1.0 - c) * pred) + (c * pii_avg)

    # get predicted and observed diversity over the chromosome
    po_arr = predicted_and_observed(cum_pos, msk, nt, dv, pred, scale, slide)

    # save predicted and observed diversity to file
    fpth = '/'.join(f_init.split('/')[:-1])
    if leave_one_out:
        fpth = root_dir + '/result/final_files/{fldr}_jackknife_results'.format(fldr=argv[1])
    if slide:
        fmt = fpth + '/{}.{:.2e}win_{:.1f}slide_pred_and_obs.txt'
        fout = fmt.format(chrom, scale, slide)
    else:
        fmt = fpth + '/{}.{:.2e}win_pred_and_obs.txt'
        fout = fmt.format(chrom, scale)
    np.savetxt(fout, po_arr)


if __name__ == '__main__':
    main()

