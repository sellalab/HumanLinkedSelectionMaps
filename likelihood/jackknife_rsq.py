__author__ = 'davidmurphy'


import numpy as np
from sys import argv
from classes.runstruct import ChromStruct, root_dir
from precalc.lh_inputs import load_saved, adjust_arrays
from likelihood.calc_rsq import calc_rsquared, standardwindows


init_dir = root_dir + '/result/init_files'
final_dir = root_dir + '/result/final_files'
ffmt = init_dir + '/YRI.{an}.BS1.6.CS0.0.NOT_STARTED.initial.txt'


def main():
    # if len(argv) != 2:
    #     print 'usage: jackknife_rsq <folder_name>'
    #     exit(1)

    if len(argv) != 2:
        print('usage: jackknife_rsq <anno>')
        exit(1)

    anno = argv[1]
    # # use the general purpose init file for the jackknife folder
    # f_init = root_dir + '/result/final_files/{}/jk_init.txt'.format(argv[1])
    # cst = ChromStruct(chrom='chr1', init=f_init)
    #
    f_init = ffmt.format(an=anno)
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
    bs, nu, nt, dv = [a[mnu2] for a in [bs, nu, nt, dv]]
    msk &= mnu1[:,0]

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

    # LOAD concatenated jackknife sample predictions
    # f_pred = root_dir + '/result/final_files/{}/jkpred.npy'.format(argv[1])
    fpred = final_dir + '/{an}_jackknife_results/{an}.jkpred.npy'.format(an=anno)
    pred = np.load(fpred)[mnu2]

    # use init file path as template for rsq file
    # rsq_file = root_dir + '/result/final_files/{}/jkrsq.log'.format(argv[1])
    rsq_file = final_dir + '/{an}_jackknife_results/{an}.jkrsq.log'.format(an=anno)

    # calculate r squared
    with open(rsq_file, 'w') as f:
        for sc in standardwindows():
            rsq = calc_rsquared(sc, cum_pos, msk, nt, dv, pred)
            f.write('{:.3e}\t{}\n'.format(sc, rsq))


if __name__ == '__main__':
    main()
