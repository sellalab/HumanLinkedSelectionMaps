__author__ = 'davidmurphy'


from shutil import move
import numpy as np
from sys import argv
from runinf import run_inference
from classes.runstruct import ChromStruct, root_dir, default_fitness_effects
from simulations import cluster_sim_config as scfg


def run_test(init, idx, min_bsx, min_b):
    # bdfe = default_fitness_effects
    # meth = 'Nelder-Mead'
    # b_anno = 'primate_cons95_Segments'
    # bdir = 'std_split_pts'
    # tkn = 'pr95.cleanrun'
    udel = 1e-08
    alph = 0.25
    # min_bs = None

    # build dict of arguments needed for cst instance
    chrom = 'chr1'
    # bs1_args = dict(bs_annos=(b_anno,), cs_annos=(),
    #                 bdfe=(bdfe,), cdfe=(), bdir=bdir,
    #                 tkn=tkn, chrom=chrom, methods=(meth,))
    # # initialize chromstruct and set simulation stats
    # cst = ChromStruct(**bs1_args)
    cst = ChromStruct(chrom=chrom, init=init)
    # cst.stat.meandiv = scfg.mean_div
    # cst.stat.indv = 216

    # set bounds for CLH function
    cst.fixed.min_bs = min_b
    cst.fixed.min_bsx = min_bsx

    # build dict of arguments needed pgm instance
    # TODO: FIX THIS PARAMETERIZATION
    nparam = cst.bsgrid if cst.bsgrid else cst.csgrid
    pgm_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=nparam,
                    itau=scfg.init_tau, ftau=scfg.true_tau,
                    udel=udel, alpha=alph)

    # initialize paramgenerator to get random input params
    pgm = scfg.ParamGenerator(**pgm_args)

    # get initial params using index given
    subset = 20, 25
    scales = (0.05, 0.5, 5.0)

    # # TODO: integrate this setting if we end up using it
    # params_1 = pgm.random_set(subset, scales)[idx]
    # lsprm = params_1[:-1]
    # tau = params_1[-1]
    # cth_init = 0.01
    # params_2 = np.concatenate((lsprm, [cth_init], [tau]))
    # cst.params = params_2

    cst.params = pgm.random_set(subset, scales)[idx]
    cst.stat.init_params = cst.params

    # get neutral sim file token and create a run ID string
    run_id = 'iprm_{:02}'.format(idx)
    bth_str = 'bth_{:.3f}'.format(min_bsx)
    if min_b is not None:
        mb_str = '_minb_{:.2e}'.format(min_b)
        bth_str += mb_str

    # rename text file to include more information about run parameters
    f = '{rslt}/init_files/{neut}.{lbl}.{i}.{n}.{timestamp}.initial.txt'
    cst.txt_file = f.format(i=run_id, n=bth_str, **cst.dict)

    # run the inference
    run_inference(cst, parallel=True)  # parallel not working on personal computer for some reason...
    # run_inference(cst, parallel=False)


def main():
    # path to initialization file for the main CADD 6% results
    init = root_dir + '/result/init_files/YRI.cadd94_gmask_v1.6_without_bstat.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    # minimum b-value to apply *within optimization* (0=None)
    min_bsx = 0
    # minimum b-value to apply to precalc maps
    min_b = -3.78010917034472
    # loop through all 15 indices and run inference for each set of initial conditions:
    for idx in range(15):
        run_test(init, idx, min_bsx, min_b)

    # if root_dir.startswith('/Users/MURPHYD/'):
    #     # initialization file (contains params, points to files, etc)
    #     init = root_dir + '/result/init_files/YRI.fish_cons94_new.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    #     # index of starting parameters to use (0-14)
    #     idx = 0
    #     # minimum b-value to apply *within optimization* (0=None)
    #     min_bsx = 0
    #     # minimum b-value to apply to precalc maps
    #     min_b = -3.78010917034472
    #
    # # take arguments from the command line (used in cluster)
    # else:
    #     if len(argv) != 5:
    #         print('usage: cluster_runinf <init> <idx> <min_bsx> <min_b')
    #         exit(1)
    #     init = argv[1]
    #     idx = int(argv[2])
    #     min_bsx = eval(argv[3])
    #     min_b = eval(argv[4])
    #     if min_b is not None:
    #         min_b = np.log(min_b)
    # run_test(init, idx, min_bsx, min_b)


if __name__ == '__main__':
    main()
