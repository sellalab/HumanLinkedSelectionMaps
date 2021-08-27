__author__ = 'davidmurphy'

import os, re
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import simulations.cluster_sim_config as scfg
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import adjust_arrays, load_saved
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir
dfe = default_fitness_effects


def make_plots(folder_name):
    bth = re.search('bth(0\.\d+)', folder_name).group(1)
    fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
    # empty list for results files
    rlst = []
    # collect results files and create ChromStructs
    for f in os.listdir(fdir):
        if f.endswith('.txt'):
            fname = fdir + f
            rlst.append(ChromStruct(chrom='chr1', init=fname))

    # get an indexed list of all likelihoods for run list
    ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
    # get indices of the top 3 best LH runs (after sorting on LH)
    top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

    # create the plot
    plt.figure(figsize=(15, 5))
    plt.suptitle('B threshold = '+bth)
    plt.subplots_adjust(left=0.08, wspace=1.5, right=0.9)
    plt.subplot(1, 12, (1, 4))
    if 'sample' in fdir:
        plt.title('background selection sampled n=2')
    else:
        plt.title('background selection deterministic')

    x = np.arange(6)
    n = 2
    w = 0.8 / n
    s = -w / 2.0

    # get the true params for the run and plot as the first variable
    pgn_args = dict(nba=rlst[0].bnum, nca=rlst[0].cnum, nprm=len(dfe),
                    itau=scfg.init_tau, ftau=scfg.true_tau,
                    udel=scfg.mean_u, alpha=0.25)
    # initialize paramgenerator to get random input params for the simulation
    pgn = scfg.ParamGenerator(**pgn_args)
    true_params = pgn.random(tau=pgn_args['ftau'])[0]

    tprm = 10 ** true_params[:-1]
    # tprm = 10 ** rlst[0].stat.true_params[:-1]
    tpmf = tprm * 1e8 * 0.9
    plt.bar(x + s, tpmf, w, color='k', label='true')
    s += w

    # plot results from the top3 indices of runs (averaged)
    pmfs = []
    for i in top3:
        r = rlst[i]
        pmfs.append(r.uvec[0])
    pmfs = np.array(pmfs)
    mpmf = np.average(pmfs, axis=0) * 1e8

    plt.bar(x + s, mpmf, w, label='best')
    plt.ylabel(r'$\mu_{del} \cdot 10^8$')
    plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
    plt.xlabel('deleterious fitness effect')
    plt.legend()

    # plot udel
    ud = sum(mpmf)  # udel
    udt = sum(tpmf) * 0.9 # true udel
    sp3 = plt.subplot(1,12,5)
    sp3.yaxis.tick_right()
    sp3.yaxis.set_label_position('right')
    plt.xlabel(r'$\mu_{del} \cdot 10^8$')
    plt.bar(0, udt, 1, color='k')
    plt.bar(1, ud, 1)
    plt.xticks([])

    # plot pi0
    plt.subplot(1,12,6)
    plt.xlabel(r'$\pi_0$')
    pi0 = rlst[top3[0]].params[-1]
    tpi0 = true_params[-1]
    plt.bar(0, tpi0, 1, color='k')
    plt.bar(1, pi0, 1)
    plt.xticks([])
    plt.yticks([])

    # plot clh_best - clh_true
    sp4 = plt.subplot(1,12,7)
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    sp4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel(r'$\Delta CLH$')
    # plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
    # clh, tclh = rlst[top3[0]].stat.best_lh, rlst[top3[0]].stat.true_clh
    clh = rlst[top3[0]].stat.best_lh
    tclh = 797.5930353871922
    diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
    col = 'green' if diff > 0 else 'red'
    plt.bar(-0.4, diff, 0.4, color=col)
    plt.ylim(0.99 * diff, 1.01 * diff)
    plt.xticks([])

    sp4 = plt.subplot(1,12, (8,12))
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    # sp4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel(r'$R^2$')
    rlog = [fdir+f for f in os.listdir(fdir) if f.endswith('.log')][0]
    rvals = np.loadtxt(rlog)
    xvals = range(3)
    plt.bar(xvals, rvals, 0.8)
    # plt.ylim(0.99 * diff, 1.01 * diff)
    plt.xlabel('window size')
    plt.xticks(xvals, [r'$10^{%.1f}$' % x for x in [4,5,6]])
    fsave = fdir + 'summary.png'
    plt.savefig(fsave, dpi=256)
    plt.close()
    # plt.show()


def compare_runs(folder_list, label_list, colors, markers, save_name):
    # create lists for data to plot
    dfe_list = []
    pi0_list = []
    clh_list = []
    rsq_list = []
    # set initial tau manually
    init_tau = 64.90887833255238

    for fldr in folder_list:
        # find composite file and rsq file in folder (should be 1 of each)
        fdir = root_dir + '/result/final_files/{}/'.format(fldr)
        flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
        rlist = [f for f in os.listdir(fdir) if f.endswith('.log')]
        assert len(flist) == len(rlist) == 1

        # create paths to rsq and final init file
        f_init = fdir + flist[0]
        r_log = fdir + rlist[0]

        # get the rsq values
        rsq = np.loadtxt(r_log)
        rsq_list.append(rsq)

        # initalize ChromStruct with init file for final params
        cst = ChromStruct('chr1', init=f_init)

        # get dfe and pi0 from composite run
        mpmf = cst.uvec[0] * 1e8
        mpi0 = cst.params[-1] / init_tau
        dfe_list.append(mpmf)
        pi0_list.append(mpi0)

        # get best CLH
        clh_list.append(cst.stat.best_lh)


    # create the plot
    plt.figure(figsize=(12, 5))
    plt.suptitle('comparing B thresholds')
    plt.subplots_adjust(left=0.08, wspace=1.5, right=0.9)
    plt.subplot(1,11, (1, 4))
    plt.title('distribution of fitness effects')

    x = np.arange(6)
    n = len(folder_list)
    w = 0.8 / n
    s = -w * n / 2.0

    # plot results from the top3 indices of runs (averaged)
    for (i, df) in enumerate(dfe_list):
        lbl = label_list[i]
        if len(df) < 6:
            nz = 6 - len(df)
            zpad = [0] * nz
            df = np.concatenate((zpad, df))
        plt.bar(x + s, df, w, label=lbl, color=colors[i], align='edge')
        s += w
    plt.ylabel(r'$\mu_{del} \cdot 10^8$')
    plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
    plt.xlabel('deleterious fitness effect')
    plt.legend(prop=dict(size=9))

    # plot udel
    sp3 = plt.subplot(1,11,5)
    sp3.yaxis.tick_right()
    sp3.yaxis.set_label_position('right')
    plt.xlabel(r'$\mu_{del} \cdot 10^8$')
    s = -w / 2.0
    for (i, df) in enumerate(dfe_list):
        udel = sum(df)
        plt.bar(0+s, udel, w, color=colors[i])
        s += w
    plt.xticks([])

    # plot pi0
    sp4 = plt.subplot(1,11,6)
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    plt.xlabel(r'$\frac{\pi}{\pi_0}$')
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=colors[i])
        s += w
    plt.xticks([])
    # plt.yticks([])

    # plot CLH
    sp4 = plt.subplot(1,11,7)
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    sp4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('-CLLH')
    s = -w * n / 2.0
    clh_list = [cl/1e5 for cl in clh_list]
    for (i, cl) in enumerate(clh_list):
        plt.bar(0+s, cl, w, color=colors[i])
        s += w
    plt.ylim(min(clh_list), max(clh_list))
    # plt.ylim(0.999*min(clh_list), 1.001*max(clh_list))
    plt.xticks([])

    sp4 = plt.subplot(1,11,(8,11))
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    plt.title('variance explained')
    plt.ylabel(r'$R^2$')
    # get window sizes from the first file
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    # s = -w * n / 2.0
    for (i, rs) in enumerate(rsq_list):
        rsq = rs[:,1]
        plt.plot(ws[wmsk], rsq[wmsk], color=colors[i], marker=markers[i],
                 linewidth=0, alpha=0.95, label=label_list[i])
        # plt.bar(x+s, rsq, w, color=colors[i], align='edge')
        # s += w
    plt.xlabel('window size')
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4,5,6]])
    plt.legend(prop=dict(size=9))

    # save the plot
    fsave = root_dir + '/result/final_files/{}.png'.format(save_name)
    # plt.savefig(fsave, dpi=256)
    # plt.close()
    plt.show()


def main():
    # set colors for bar plots
    # clist = 'hotpink deeppink mediumvioletred mediumseagreen ' \
    #         'darkgreen aquamarine darkturquoise darkcyan'.split()
    # flist = []
    # llist = []
    # for bth in [0, 650]:
    #     folder_name = 'datarun_{:03}'.format(bth)
    #     label_name = r'$B \geq {}$'.format(float(bth) / 1000)
    #     flist.append(folder_name)
    #     llist.append(label_name)
    #
    # flist += ['std_run_feb2019', 'cut_dfe1', 'cut_dfe2', 'nff2', 'nff2_650',
    #           'nff2_bsmin010']
    # llist += [r'$precalc\ B \geq 0.01$', r'$t \geq 10^{-3.5}$',
    #           r'$t \geq 10^{-3}$', r'$B \geq 0\ old\ data$',
    #           r'$B \geq 0.65\ old\ data$', r'$precalc\ B \geq 0.01\ old\ data$']
    # mlist = 'o o o ^ ^ s s s'.split()
    # flist = 'datarun_000 nff2 datarun_650 nff2_650 std_run_feb2019 ' \
    #         'nff2_bsmin010'.split()
    # llist = [r'$B \geq 0$', r'$B \geq 0\ (old)$', r'$B \geq 0.65$',
    #          r'$B \geq 0.65\ (old)$', r'$precalc\ B \geq 0.01$',
    #          r'$precalc\ B \geq 0.01\ (old)$']
    # clist = 'deeppink mediumvioletred darkturquoise darkcyan goldenrod ' \
    #         'chocolate'.split()
    # mlist = 'o o s s ^ ^'.split()
    # sname = 'compare_oldnew'
    # flist = 'datarun_000 datarun_650 std_run_feb2019 cut_dfe1 cut_dfe2'.split()
    # llist = [r'$B \geq 0$', r'$B \geq 0.65$', r'$precalc\ B \geq 0.01$',
    #          r'$t \geq 10^{-3.5}$', r'$t \geq 10^{-3}$']
    # clist = 'hotpink deeppink mediumvioletred mediumseagreen darkgreen'.split()
    # mlist = 'o o o ^ ^'.split()
    # sname = 'compare_truncated_dfe'
    flist = 'std_run_feb2019 nff2_bsmin010 rmneut_bsmin010'.split()
    llist = [r'$precalc\ B \geq 0.01\ remove\ cons\ neut$',
             r'$precalc\ B \geq 0.01\ keep\ cons$',
             r'$precalc\ B \geq 0.01\ remove\ neut\ cons$']
    clist = 'hotpink mediumseagreen goldenrod'.split()
    mlist = 'o s ^'.split()
    sname = 'compare_cons_filters'
    compare_runs(flist, llist, clist, mlist, sname)


main()