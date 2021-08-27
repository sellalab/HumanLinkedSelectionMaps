import os
import numpy as np
from itertools import izip
from collections import defaultdict, namedtuple
from sys import argv, stderr, stdout
from shutil import move
from likelihood.runinf import run_inference, evaluate
from classes.runstruct import ChromStruct, RunStruct, root_dir, \
    default_fitness_effects
from simulations import cluster_sim_config as scfg
dfe = default_fitness_effects


__author__ = 'davidmurphy'


def results_table(fdir):
    """
    Compare inferred params vs. true params from simulated data using
    different optimizers under varying sets of initial conditions against
    """
    # add "/" to final file directory if needed for fdir + f concatenation
    fdir = fdir if fdir.endswith('/') else fdir + '/'

    # collect the optimization run results stored in the directory fdir
    files = [fdir + f for f in os.listdir(fdir) if 'iprm_' in f]

    # get initial conditions indices from run file names
    idx = [int(f.split('iprm_')[-1][:2]) for f in files]

    # initialize runstructs from each file
    rs = [RunStruct(init=f) for f in files]

    # create sparse array for results using nested dict[idx][method] = result)
    res_tab = defaultdict(dict)

    # create namedtuple to organize key results from run
    Result = namedtuple('Result', 'prm clh runtime fcalls status isnan')

    # store results from each file in the table
    for (i, r) in izip(idx, rs):
        # inferred params
        prm = r.params
        # CLH at initial and inferred params
        # clh = (r.stat.initial_lh, r.op1.fun)
        clh = r.op1.fun
        # run time
        rt = r.stat.total_time
        # function calls
        fcalls = r.stat.function_calls
        # optimization exit status
        stat = int(r.op1.success)
        # nan flag for final CLH value
        isnan = np.isnan(r.op1.fun)
        # get the method for the present result
        meth = r.methods[0]

        # store results in the appropriate cells of res_tab
        res_tab[i][meth] = Result(prm, clh, rt, fcalls, stat, isnan)

    return res_tab


def run_test_1():
    """
    Usage: set initial parameters manually and run parallelized optimization on the
    cluster.
    """
    if root_dir.startswith('/Users/davidmurphy'):
        # iprm_idx = 1
        meth = ('Nelder-Mead',)
        dst = 'rnd_0'
        smp = 2
        idx = 0
        itr = 0

    else:
        if len(argv) != 5:
            # print 'usage: convergence_tests <meth> <dist> <fact> <rand>'
            # print 'usage: convergence_tests <meth> <dist> <idx>'
            # print 'usage: convergence_tests <meth> <dist>'
            print 'usage: convergence_tests <meth> <dist> <sample> <idx>'
            exit(1)

        # iprm_idx = int(argv[1])
        meth = tuple(argv[1].split(','))
        dst = argv[2]
        smp = int(argv[3])
        idx = int(argv[4])
        itr = 0

        # stderr.write(meth[0]+'\n')
        # stdout.flush()
    # build dict of arguments needed for cst instance
    cst_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(scfg.c_anno,),
                    methods=meth, bdfe=(dfe,), cdfe=(dfe,), bdir=scfg.bdir,
                    cdir=scfg.cdir, tkn=scfg.tkn, chrom='chr1')

    # initialize chromstruct and set simulation stats
    cst = ChromStruct(**cst_args)
    cst.stat.meandiv = scfg.mean_div
    cst.stat.indv = smp

    # # create new chrom struct from config params
    # cst = ChromStruct(chrom='chr1', bdir=scfg.bdr, tkn=scfg.tkn,
    #                   methods=meth)

    # # from data_processing.functions import swap_root
    # f = idict[(dst, rndm)]
    # # swap_root(f)
    # stderr.write('use init file={}\n'.format(f))
    # stdout.flush()
    # old_cst = ChromStruct(chrom='chr1', init=f)

    # set bounds for CLH function
    cst.fixed.min_bs = scfg.min_bs
    cst.fixed.min_red = scfg.min_red

    # # SET INITIAL PARAMS FROM AVERAGE 5 BEST NM PARAMS
    # rdir = '{}/result/final_files/opttests/bs_cs'.format(root_dir)
    # if not os.path.isdir(rdir):
    #     os.mkdir(rdir)
    # ridx = int(dst.split('_')[-1])
    # td = '{}/test_{}'.format(direct, tidx+5)
    # fs = ['{}/{}'.format(td, f) for f in os.listdir(td) if 'YRI' in f]
    # cs = [ChromStruct('chr1', init=f) for f in fs]
    # sort_cs = sorted(cs, key=lambda c: c.stat.best_lh)
    # p_init = np.average(np.array([c.params for c in sort_cs[:3]]), axis=0)

    # set initial parameters
    # mean_u, itau = scfg.mean_u, scfg.init_tau
    # p_init = scfg.uniform_params(mean_u, fact, itau, rndm)
    # p_init = np.array(scfg.iprm[idx_idx])

    # build dict of arguments needed pgm instance
    pgm_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=len(dfe),
                    itau=scfg.init_tau, ftau=scfg.true_tau, udel=scfg.mean_u,
                    alpha=0.25)

    # initialize paramgenerator to get random input params for the simulation
    pgm = scfg.ParamGenerator(**pgm_args)

    # get true param row index from dst label
    ridx = int(dst.split('_')[-1])
    tprm = pgm.random(tau=pgm_args['ftau'])[ridx]

    # select random initial at idx
    subset = 20, 25
    scales = (0.05, 0.5, 5.0)
    iprm = pgm.random_set(subset, scales)[idx]

    cst.params = iprm
    cst.init_params = cst.params

    # msg = 'init_params = {}'.format(str(cst.params))
    # stderr.write(msg+'\n')
    # stdout.flush()

    # set nt files to sim names
    # newtkn, tprm = scfg.true_params(dst)
    # newtkn = newtkn.replace(dst, dst+'_noisy')
    # samp = 100
    # newtkn = newtkn.replace(dst, dst+'_n={}'.format(samp))

    # get neutral sim file token
    newtkn = 'nuconst_{}_smpl_{}_itr_{:03}'.format(dst, smp, itr)
    cst.files.fnt = cst.nt_simfiles(newtkn)
    cst.true_params = tprm

    # create new text file name based on index and method
    f = '{rslt}/init_files/{neut}.{lbl}.{i}.{m}.{n}.{timestamp}.initial.txt'

    # run ID label
    # run_id = 'iprm_{:02}'.format(iprm_idx)
    run_id = 'mean_NM'
    # run_id = 'init_ufm_udel'
    # if rndm:
    #     run_id += '_rndm'

    # methods label
    # methods = '_'.join(meth)
    # methods = 'NM_' + meth[0]
    methods = meth[0]

    # new text file
    cst.txt_file = f.format(i=run_id, m=methods, n=newtkn, **cst.dict)

    # set bounds manually
    if meth[0] == 'trust-constr':
        cst.op1.bounds = cst.fixed.min_max_bounds
        cst.op1.bounds[-1] = (cst.params[-1]*0.1, cst.params[-1]*10)
        cst.op1.jac = True
        cst.op1.hess = '2-point'
        cst.vars.calc_deriv = True

    # from scipy.optimize import BFGS, HessianUpdateStrategy, SR1
    # cst.op1.hess = BFGS(exception_strategy='skip_update')
    # cst.op1.hess.exception_strategy = 'damp_update'

    run_inference(cst, parallel=True)

    # # timestamp init file and save (timestamp format=YYMMDDHHMMSS)
    # cst.record_timestamp()
    # # cst.save()
    #
    # # switch off the 'complete' flag before entering optimization
    # cst.vars.complete = False
    #
    # # loop through optimization methods
    # for optimizer_id in cst.optimizers:
    #     # select optimizer
    #     optimizer = cst[optimizer_id]
    #
    #     # # DEBUG
    #     # optimizer.options['maxfev'] = 1
    #     # optimizer.options['maxiter'] = 1
    #
    #     # set derivative flag based on the method
    #     if optimizer.method in scfg.grad_methods:
    #         cst.vars.calc_deriv = True
    #         optimizer.jac = True
    #     else:
    #         cst.vars.calc_deriv = False
    #         optimizer.jac = None
    #
    #     # run the inference using the selected optimizer
    #     msg = 'starting {} optimization\n'.format(optimizer.method)
    #     stderr.write(msg)
    #     stdout.flush()
    #
    #     run_inference(cst, optimizer, parallel=True)
    #
    #     # write the final output to the optimization log
    #     if len(cst.cached_params) > 0:
    #         cst.log_params()
    #
    #     # mark the run as complete and calculate final stats
    #     cst.vars.complete = True
    #     cst.stat.calc_stats(cst)
    #
    #     # save the final file
    #     cst.save(txt_file=cst.final_file)
    #
    #     msg = 'finished {} optimization\n'.format(optimizer.method)
    #     stderr.write(msg)
    #     stdout.flush()

    # # write the final output to the optimization log
    # if len(cst.cached_params) > 0:
    #     cst.log_params()
    #
    # # mark the run as complete and calculate final stats
    # cst.vars.complete = True
    # cst.stat.calc_stats(cst)
    #
    # # save the final file
    # cst.save(txt_file=cst.final_file)


def run_test_2():
    """
    Usage: set initial parameters manually and run parallelized optimization on the
    cluster.
    """
    if root_dir.startswith('/Users/davidmurphy'):
        init = root_dir + '/result/init_files/YRI.pr95.cleanrun.ds.01.BS1.6.CS0.0.NOT_STARTED.initial.txt'
        meth = ('Nelder-Mead',)
        dst = 'rnd_0'
        smp = 2
        idx = 0
        rnd = True
        fudel = 1
        alpha = 0.25
        bthresh = 0.4

    else:
        if len(argv) != 10:
            print 'usage: convergence_tests <init> <meth> <dist> <sample> ' \
                  '<idx> <rnd> <fudel> <alpha> <bthresh>'
            exit(1)

        init = argv[1]
        meth = tuple(argv[2].split(','))
        dst = argv[3]
        smp = int(argv[4])
        idx = int(argv[5])
        rnd = eval(argv[6])
        fudel = float(argv[7])
        alpha = float(argv[8])
        bthresh = float(argv[9])

    # iteration tag
    itr = 0

    # build dict of arguments needed for cst instance
    # bs1cs1_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(scfg.c_anno,),
    #                    bdfe=(dfe,), cdfe=(dfe,), bdir=scfg.bdir,
    #                    cdir=scfg.cdir, tkn=scfg.tkn, chrom=chrom, methods=meth)
    # cs1_args = dict(bs_annos=(), cs_annos=(scfg.c_anno,),
    #                 bdfe=(), cdfe=(dfe,), cdir=scfg.cdir,
    #                 tkn=scfg.tkn, chrom=chrom, methods=meth)
    # bs1_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(),
    #                 bdfe=(dfe,), cdfe=(), bdir=scfg.bdir,
    #                 tkn=scfg.tkn, chrom=chrom, methods=meth)
    #
    # # set the arg dict to be used
    # if stype == 'bs1':
    #     use_args = bs1_args
    # elif stype == 'cs1':
    #     use_args = cs1_args
    # else:
    #     assert stype == 'bs1cs1'
    #     use_args = bs1cs1_args

    # initialize chromstruct and set simulation stats
    cst = ChromStruct(chrom='chr1', init=init)
    assert cst.methods == meth
    cst.stat.meandiv = scfg.mean_div
    cst.stat.indv = smp

    # set bounds for CLH function
    cst.fixed.min_bs = scfg.min_bs
    cst.fixed.min_red = scfg.min_red
    cst.fixed.min_bsx = bthresh

    # set constant mutation rate flag in variables to true
    cst.vars.mu_const = True

    # build dict of arguments needed pgn instance
    pgn_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=len(dfe),
                    itau=scfg.init_tau, ftau=scfg.true_tau,
                    udel=scfg.mean_u, alpha=alpha)

    # initialize paramgenerator to get random input params for the simulation
    pgn = scfg.ParamGenerator(**pgn_args)

    # get random param row index from dst label
    if dst.startswith('rnd'):
        ridx = int(dst.split('_')[-1])
        cst.stat.true_params = pgn.random(tau=pgn_args['ftau'])[ridx]
    elif dst.startswith('ufm'):
        cst.stat.true_params = pgn.uniform(pattern=dst, tau=pgn_args['ftau'])
    else:
        raise NameError('{} distribution does not exist'.format(dst))

    # get initial params using index given
    subset = 20, 25
    scales = (0.05, 0.5, 5.0)
    cst.params = pgn.random_set(subset, scales)[idx]
    cst.stat.init_params = cst.params

    # get neutral sim file token and create a run ID string
    run_id = 'iprm_{:02}'.format(idx)
    # newtkn = 'nuconst_{}_smpl_{}_itr_{:03}'.format(dst, smp, itr)
    stkn = 'nuconst_{}_smpl_{}_itr_{:03}_alpha_{:.1e}_fudel_{:.1e}'
    newtkn = stkn.format(dst, smp, itr, alpha, fudel)
    if not rnd:
        newtkn += '_detrm'
    cst.files.fnt = cst.nt_simfiles(newtkn)

    # rename text file to include more information about run parameters
    f = '{rslt}/init_files/{neut}.{lbl}.{i}.{m}.{n}.{timestamp}.initial.txt'
    cst.txt_file = f.format(i=run_id, m=meth[0], n=newtkn, **cst.dict)

    # set bounds manually if gradient method is used
    if meth[0] == 'trust-constr':
        cst.op1.bounds = cst.fixed.min_max_bounds
        cst.op1.bounds[-1] = (cst.params[-1]*0.1, cst.params[-1]*10)
        cst.op1.jac = True
        cst.op1.hess = '2-point'
        cst.vars.calc_deriv = True
    if meth[0] in ('CG', 'BFGS'):
        cst.op1.jac = True
        cst.vars.calc_deriv = True

    # run the inference
    run_inference(cst, parallel=True)


def run_test_3():
    """
    Usage: set initial parameters manually and run parallelized optimization on the
    cluster.
    """
    if root_dir.startswith('/Users/davidmurphy'):
        init = root_dir + '/result/init_files/YRI.pr95.cleanrun.ds.01.BS1.6.CS0.0.NOT_STARTED.initial.txt'
        idx = 0
        bthresh = 0.4

    else:
        if len(argv) != 4:
            print 'usage: convergence_tests <init> <idx> <bthresh>'
            exit(1)
        init = argv[1]
        idx = int(argv[2])
        bthresh = float(argv[3])

    # initialize chromstruct and set simulation stats
    cst = ChromStruct(chrom='chr1', init=init)
    # check assumptions for the run
    assert cst.methods == ('Nelder-Mead',)
    # assert cst.stat.indv == 2
    cst.stat.meandiv = scfg.mean_div
    cst.stat.indv = 2

    # set bounds for CLH function
    cst.fixed.min_bsx = bthresh

    # set constant mutation rate flag in variables to true
    cst.vars.mu_const = True

    # build dict of arguments needed pgn instance
    alpha = 0.25
    pgn_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=len(dfe),
                    itau=scfg.init_tau, ftau=scfg.true_tau,
                    udel=scfg.mean_u, alpha=alpha)

    # initialize paramgenerator to get random input params for the simulation
    pgn = scfg.ParamGenerator(**pgn_args)

    # get initial params using index given
    subset = 20, 25
    scales = (0.05, 0.5, 5.0)
    cst.params = pgn.random_set(subset, scales)[idx]
    cst.stat.init_params = cst.params

    # get neutral sim file token and create a run ID string
    run_id = 'iprm_{:02}'.format(idx)
    # # newtkn = 'nuconst_{}_smpl_{}_itr_{:03}'.format(dst, smp, itr)
    # stkn = 'nuconst_{}_smpl_{}_itr_{:03}_alpha_{:.1e}_fudel_{:.1e}'
    # newtkn = stkn.format(dst, smp, itr, alpha, fudel)
    # if not rnd:
    #     newtkn += '_detrm'
    # cst.files.fnt = cst.nt_simfiles(newtkn)
    print 'USING NEUTRAL POLYMORPHISM FROM FILE: {}'.format(cst.nt_files)

    # # rename text file to include more information about run parameters
    f = '{rslt}/init_files/{neut}.{lbl}.{i}.{timestamp}.initial.txt'
    cst.txt_file = f.format(i=run_id, **cst.dict)

    # run the inference
    run_inference(cst, parallel=True)

def main():
    run_test_3()


if __name__ == '__main__':
    main()
