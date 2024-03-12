import os
import numpy as np
from sys import argv
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import load_saved, adjust_arrays
from likelihood.calc_rsq import calc_rsquared, standardwindows
from classes.runstruct import root_dir, cst_from_fldr, human_autosomes


init_dir = root_dir + '/result/init_files'
final_dir = root_dir + '/result/final_files'
ffmt = init_dir + '/YRI.{an}.BS1.6.CS0.0.NOT_STARTED.initial.txt'


def dropped_chrom_prediction():
    """
    Create a prediction vector by combining chromosomes, where predictions
    for each chromosome are based on the inference run without data from
    that chromosome.
    """
    fldr = 'cadd94_gmask_v1.6_without_bstat'
    mean_div = cst_from_fldr(fldr).stat.meandiv
    pred_list = []
    pred_const_list = []
    # for each chromosome, load the results with that chromosome dropped
    for ch in human_autosomes:
        # load results with the current chrom dropped
        fldr = 'drop_{}'.format(ch)
        cst = cst_from_fldr(fldr)
        # load the maps for the chrom that was excluded
        cst.vars.drop_chrom = None
        sg, bs, cs, nu, nt, dv, pl = load_saved(cst, chroms=[ch])
        bs, cs, nu = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)[:3]
        # get constant mutation rate for sorted maps
        nu_const = np.full(len(nu), mean_div)
        # generate predictions for params inferred with the chrom excluded
        pred = predicted_pi(cst.params, cst, nu, bs, cs)
        pred_const = predicted_pi(cst.params, cst, nu_const, bs, cs)
        # add predictions to list from other chroms
        pred_list.append(pred)
        pred_const_list.append(pred_const)

    # concatenate all of the predictions for the 22 autosomes and save
    if not os.path.isdir(final_dir + '/drop_chrom_results'):
        os.mkdir(final_dir + '/drop_chrom_results')
    f_pred = final_dir + '/drop_chrom_results/drop_chrom_pred.npy'
    f_pred_const = final_dir + '/drop_chrom_results/drop_chrom_pred_nuconst.npy'
    np.save(f_pred, np.concatenate(pred_list))
    np.save(f_pred_const, np.concatenate(pred_const_list))


def dropped_chrom_rsq():
    """
    Calculate R^2 using droppped chroms prediction map
    """
    fldr = 'cadd94_gmask_v1.6_without_bstat'
    cst = cst_from_fldr(fldr)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    cutoff = 0.25
    mnu1 = (nu <= cutoff)

    # mask and rescale maps
    bs, cs, nu, nt, dv, pl, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    mnu2 = (nu <= cutoff)
    bs, nu, nt, dv = [a[mnu2] for a in [bs, nu, nt, dv]]
    msk &= mnu1[:, 0]

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

    # LOAD dropped chrom predictions array
    f_pred = final_dir + '/drop_chrom_results/drop_chrom_pred.npy'
    pred = np.load(f_pred)[mnu2]

    # use init file path as template for rsq file
    rsq_file = final_dir + '/drop_chrom_results/drop_chrom_rsq.log'

    # calculate r squared
    with open(rsq_file, 'w') as f:
        for sc in standardwindows():
            rsq = calc_rsquared(sc, cum_pos, msk, nt, dv, pred)
            f.write('{:.3e}\t{}\n'.format(sc, rsq))


def main():
    f_pred = final_dir + '/drop_chrom_results/drop_chrom_pred.npy'
    f_pred_const = final_dir + '/drop_chrom_results/drop_chrom_pred_nuconst.npy'
    if (not os.path.isfile(f_pred)) or (not os.path.isfile(f_pred_const)):
        dropped_chrom_prediction()
    rsq_file = final_dir + '/drop_chrom_results/drop_chrom_rsq.log'
    if not os.path.isfile(rsq_file):
        dropped_chrom_rsq()


if __name__ == '__main__':
    main()
