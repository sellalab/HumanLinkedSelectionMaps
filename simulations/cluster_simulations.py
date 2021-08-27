__author__ = 'davidmurphy'


from sys import argv
import cluster_sim_config as scfg
from neutpoly_sim import run_simulation
from classes.runstruct import ChromStruct, root_dir, default_fitness_effects
dfe = default_fitness_effects


def main():
    if root_dir.startswith('/Users/davidmurphy'):
        init = root_dir + '/result/init_files/YRI.pr95.cleanrun.ds.00.BS1.6.CS0.0.NOT_STARTED.initial.txt'
        ch = 'chr22'
        dst = 'rnd_0'
        smp = 2
        rnd = True
        fudel = 1.0
        alpha = 0.25
        bthresh = 0.4
    else:
        if len(argv) != 9:
            print 'usage: cluster_simulations <init> <ch> <dist> <sample> <random> <f_udel> <alpha> <bthresh>'
            exit(1)
        init = argv[1]
        ch = argv[2]
        dst = argv[3]
        smp = int(argv[4])
        rnd = eval(argv[5])
        fudel = float(argv[6])
        alpha = float(argv[7])
        bthresh = float(argv[8])

    # iteration tag
    itr = 0

    # # build dict of arguments needed for cst instance
    # bs1cs1_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(scfg.c_anno,),
    #                    bdfe=(dfe,), cdfe=(dfe,), bdir=scfg.bdir,
    #                    cdir=scfg.cdir, tkn=scfg.tkn, chrom=ch)
    # cs1_args = dict(bs_annos=(), cs_annos=(scfg.c_anno,),
    #                 bdfe=(), cdfe=(dfe,), cdir=scfg.cdir,
    #                 tkn=scfg.tkn, chrom=ch)
    # bs1_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(),
    #                 bdfe=(dfe,), cdfe=(), bdir=scfg.bdir,
    #                 tkn=scfg.tkn, chrom=ch)
    #
    # use_args = bs1_args


    # initialize chromstruct and set simulation stats
    # cst = ChromStruct(**use_args)
    cst = ChromStruct(chrom=ch, init=init)
    cst.stat.meandiv = scfg.mean_div
    cst.stat.indv = smp
    cst.fixed.min_bsx = bthresh

    # set constant mutation rate flag in variables to true
    cst.vars.mu_const = True

    # set alpha and udel params
    udel = scfg.mean_u * fudel

    # build dict of arguments needed pgm instance
    pgm_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=len(dfe),
                    itau=scfg.init_tau, ftau=scfg.true_tau, udel=udel,
                    alpha=alpha)

    # initialize paramgenerator to get random input params for the simulation
    pgm = scfg.ParamGenerator(**pgm_args)

    # get random param row index from dst label
    if dst.startswith('rnd'):
        ridx = int(dst.split('_')[-1])
        tprm = pgm.random(tau=pgm_args['ftau'])[ridx]
    elif dst.startswith('ufm'):
        tprm = pgm.uniform(pattern=dst, tau=pgm_args['ftau'])
    else:
        raise NameError('{} distribution does not exist'.format(dst))

    # get base file name and true params, modify filename based on sample
    stkn = 'nuconst_{}_smpl_{}_itr_{:03}_alpha_{:.1e}_fudel_{:.1e}'
    newtkn = stkn.format(dst, smp, itr, alpha, fudel)
    if not rnd:
        newtkn += '_detrm'

    # run simulations
    run_simulation(cst, newtkn, smp, params=tprm, random=rnd)


if __name__ == '__main__':
    main()
