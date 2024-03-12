__author__ = 'davidmurphy'


import numpy as np
from sys import argv
from runinf import run_inference
from classes.runstruct import ChromStruct
from simulations import cluster_sim_config as scfg


def run_jackknife_test(init, idx, jidx, min_bsx=0):
    udel = 1e-08
    alph = 0.25

    # build dict of arguments needed for cst instance
    chrom = 'chr1'
    # # initialize chromstruct and set simulation stats
    cst = ChromStruct(chrom=chrom, init=init)

    # set bounds for CLH function
    cst.fixed.min_bsx = min_bsx
    cst.fixed.min_bs = np.log(0.6**7.4)

    # set use_jackknife to TRUE and set index with default 2Mb window
    cst.vars.use_jackknife = True
    cst.vars.jackknife_index = jidx
    cst.vars.jackknife_window = 2e6

    # build dict of arguments needed pgm instance
    pgm_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=cst.bsgrid,
                    itau=scfg.init_tau, ftau=scfg.true_tau,
                    udel=udel, alpha=alph)

    # initialize paramgenerator to get random input params
    pgm = scfg.ParamGenerator(**pgm_args)

    # get initial params using index given
    subset = 20, 25
    scales = (0.05, 0.5, 5.0)
    cst.params = pgm.random_set(subset, scales)[idx]
    cst.stat.init_params = cst.params

    # create iprm, b thresh and jackkinfe idx strings for file name
    run_id = 'iprm_{:02}'.format(idx)
    bth_str = 'bth_{:.3f}'.format(min_bsx)
    jck_str = 'jkidx_{:04}'.format(jidx)

    # rename text file to include more information about run parameters
    # f = '{rslt}/init_files/{neut}.{lbl}.{i}.{n}.{j}.{timestamp}.initial.txt'
    # cst.txt_file = f.format(i=run_id, n=bth_str, j=jck_str, **cst.dict)

    f = '{rslt}/init_files/{neut}.{lbl}.{i}.{j}.{timestamp}.initial.txt'
    cst.txt_file = f.format(i=run_id, j=jck_str, **cst.dict)

    # run the inference
    run_inference(cst, parallel=True)


def main():

    if len(argv) != 4:
        print('usage: cluster_runinf <init> <idx> <jackidx>')
        exit(1)
    init = argv[1]
    idx = int(argv[2])
    jidx = int(argv[3])

    run_jackknife_test(init, idx, jidx)


if __name__ == '__main__':
    main()
