__author__ = 'davidmurphy'

# import statements
import os, re
import seaborn
import numpy as np
import matplotlib as mpl
from itertools import izip
from collections import namedtuple, defaultdict
from classes.runstruct import ChromStruct, RunStruct
from simulations.cluster_sim_config import true_params, random_inits
from simulations.convergence_tests import results_table, root_dir

# from likelihood.cllh_functions import predicted_pi, serial_cllh
# from classes.runstruct import ChromStruct, root_dir
# from precalc.lh_inputs import load_saved, adjust_arrays, prepare_inputs
# from data_processing.data_tools import neutral_snps, snpcount, chromosome_length
# from data_processing.functions import calc_pi

# set high DPI on figures globally
mpl.rcParams['figure.dpi'] = 200
plt = mpl.pyplot

# struture to hold the basic variables of interest between runs
InfRes = namedtuple('InfRes', 'prm clh isnan')

# plot colors and markers to choose from
# 'OrangeRed',
# 'HotPink', 'deeppink', 'MediumVioletRed'
cols = ['orange', 'coral', 'orangered', 'deepskyblue', 'dodgerblue', 'Teal',
        'MediumAquamarine', 'DarkGreen', 'Indigo', 'MediumOrchid', 'FireBrick',
        'olive', 'rosybrown', 'navy']
marks = ['o', 'v', '^', '<', '>', 'h', '+', 'D', '8', 's', 'p', '*', 'x', 'd']

direct = '{}/result/final_files/opttests/round_5'.format(root_dir)
lh = [(773.066299395, 773.066063708),
      (813.0775224, 813.077256149),
      (733.553699545, 733.55347627),
      (785.120853133, 785.120682281),
      (793.97398652, 793.973643857),
      (785.165137126, 785.164576032),
      (729.980300284, 729.980027083),
      (770.26940153, 770.26899581),
      (713.441768521, 713.446292383)]


def index_test_results(fdir):
    """group result by simulation params and initial optimization params"""
    # get matching files in the directory
    fnames = [fdir+f for f in os.listdir(fdir) if f.endswith('final.txt')]
    # get indices of sim params and initial opt params from the filenames
    sidx = [int(f.split('_rnd_')[1].split('_smpl_')[0]) for f in fnames]
    iidx = [int(f.split('iprm_')[1][:2]) for f in fnames]

    return izip(fnames, sidx, iidx)


def plot_sorted_results(fdir, txt):
    """plot results optimization tests for each simulation"""
    # param plot colors
    # bcol = 'k darkgoldenrod goldenrod gold khaki lemonchiffon'.split()
    ccol = 'mediumblue royalblue steelblue dodgerblue deepskyblue'.split()

    # the 3 initial parameter scales used
    scls = 0.05, 0.5, 5.0

    # subdivide data by simulation and subdivide simulation data by iprm
    res = defaultdict(lambda: defaultdict(list))
    for (f, si, ii) in index_test_results(fdir):
        igrp = int(ii / 5)  # group iprms by (0,4)=0, (5,9)=1, (10,14)=2
        res[si][igrp].append(ChromStruct('chr1', init=f))

    # create 3 panels of init conditions for each simulation
    for skey in sorted(res.keys()):
        cs = []
        for ikey in res[skey]:
            # results subset
            cs += res[skey][ikey]
            # scale for init params subset
            scl = scls[ikey]

        # true params and clh for simulation
        tprm = cs[0].stat.true_params
        tclh = cs[0].stat.true_clh
        assert all(np.all(c.stat.true_params == tprm) for c in cs)
        assert all(c.stat.true_clh == tclh for c in cs)

        # inferred params and CLH for each test
        prms = np.vstack([c.params for c in cs])
        clh = np.array([c.stat.best_lh for c in cs])
        ibest = np.where(clh == clh.min())[0][0]
        # clh = np.array([tclh] + [c.stat.best_lh for c in cs])

        # number of bs & cs params in result
        nbs, ncs = cs[0].bsparams, cs[0].csparams

        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(hspace=0.3)
        # format bar width based on number of param sets including true
        # w = 0.8 / (len(prms)+1)
        w = 0.4

        # plot bs params on top
        plt.subplot(2, 8, (1, 6))
        # txt = 'BS+CS simulation {} initial param scale {}'
        # plt.suptitle(txt.format(skey, scl))
        # txt = 'BS+CS simulation {}'
        plt.suptitle(txt.format(skey))

        # start x axis offset
        x = np.arange(nbs) - 0.4
        b = (10 ** tprm[:nbs]) * 1e8
        plt.bar(x, b, width=w, align='edge', color='k', label='true')
        x += w
        b = (10 ** prms[ibest][:nbs]) * 1e8
        plt.bar(x, b, width=w, align='edge', color=ccol[2])
        # for i, p in enumerate(prms):
        #     x += w
        #     b = 10**p[:nbs]
        #     plt.bar(x, b, width=w, align='edge', color=ccol[i])
        plt.ylabel(r'$u_{del}\ \cdot 10^{-8}$')
        xlbl = [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)]
        plt.xticks(np.arange(nbs), xlbl)
        plt.xlabel('deleterious fitness effect')
        plt.legend()

        # plot cs params on bottom
        plt.subplot(2, 8, (9, 14))
        x = np.arange(ncs) - 0.4
        c = 10 ** tprm[nbs:nbs + ncs]
        plt.bar(x, c, width=w, align='edge', color='k', label='true')
        x += w
        c = 10 ** prms[ibest][nbs:nbs + ncs]
        plt.bar(x, c, width=w, align='edge', color=ccol[2])
        # for i, p in enumerate(prms):
        #     x += w
        #     c = 10 ** p[nbs:nbs+ncs]
        #     plt.bar(x, c, width=w, align='edge', color=ccol[i])
        plt.ylabel(r'$\alpha$')
        xlbl = [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)]
        plt.xticks(np.arange(ncs), xlbl)
        plt.xlabel('adaptive fitness effect')
        plt.legend()

        # plot pi0 on the right
        plt.subplot(187)
        plt.xlabel(r'$\pi_0$')
        x = -0.4
        plt.bar(x, tprm[-1], width=w, align='edge', color='k', label='true')
        x += w
        pi0 = prms[ibest][-1]
        plt.bar(x, pi0, width=w, align='edge', color=ccol[2])
        # for i, p in enumerate(prms):
        #     x += w
        #     pi0 = p[-1]
        #     plt.bar(x, pi0, width=w, align='edge', color=ccol[i])
        plt.yticks([])
        plt.xticks(color='white')

        # plot log10 of the ratio of true clh to inferred
        plt4 = plt.subplot(188)
        # plt.xlabel(r'$log_{10}(\frac{LH_{opt}}{LH_{true}})$')
        plt.xlabel(r'$LH_{opt}-LH_{true}$')

        true_clh = np.exp(-1e-05 * tclh)
        x = -0.4
        cl = np.exp(-1e-05 * clh[ibest])
        # log_ratio = np.log10(cl / true_clh)
        diff = cl - true_clh
        plt.bar(x, diff, width=w, align='edge', color=ccol[2])
        x += w
        # for i, p in enumerate(clh):
        #     cl = np.exp(-1e-05 * p)
        #     log_ratio = np.log10(cl / true_clh)
        #     plt.bar(x, log_ratio, width=w, align='edge', color=ccol[i])
        #     x += w
        # plt.yticks(color='white')
        plt4.yaxis.tick_right()
        plt.xticks(color='white')

        plt.savefig(fdir + 'fig_{}'.format(skey), dpi=256)
        plt.close()
        # plt.show()


fd = '{}/result/final_files/opttests/bs1cs1/'.format(root_dir)
tit = 'BS+CS sample=2'
txt = tit + ' simulation {}'
plot_sorted_results(fd, txt)


def some_plots():
    for t in xrange(1, 10):
        # PLOTS 1
        td = '{}/test_{}'.format(direct, t+5)
        fs = ['{}/{}'.format(td, f) for f in os.listdir(td) if 'YRI' in f]
        cs = [ChromStruct('chr1', init=f) for f in fs]
        sort_cs = sorted(cs, key=lambda c: c.stat.best_lh)
        prm = np.average(np.array([c.params for c in sort_cs[:3]]), axis=0)
        tprm = random_inits(t)
        tclh  = cs[0].stat.true_clh
        clh = np.array([c.stat.best_lh for c in sort_cs[:3]]).mean()
        # tclh, clh = lh[t-1]

        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(wspace=0.055)
        plt.subplot(1, 8, (1, 6))
        plt.title('random params {}'.format(t))
        x = np.arange(6)
        ty = 10**tprm[:-1]*1e8
        iy = 10**prm[:-1]*1e8
        plt.bar(x-0.4, ty, 0.4, label='true params', color='black')
        plt.bar(x, iy, 0.4, label='avg. top 3 NM', color='darkorange')
        plt.ylabel(r'$u_{del}\ \cdot 10^{-8}$')
        plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
        plt.xlabel('deleterious fitness effect')
        plt.legend()

        plt.subplot(187)
        plt.title('pi0')
        plt.bar(-0.4, tprm[-1], 0.4, color='black')
        plt.bar(0, prm[-1], 0.4, color='darkorange')
        plt.xticks([])
        plt.yticks([])

        sp3 = plt.subplot(188)
        sp3.yaxis.tick_right()
        sp3.yaxis.set_label_position('right')
        plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
        diff = np.exp(-1e-05 * clh)-np.exp(-1e-05 * tclh)
        col = 'green' if diff > 0 else 'red'
        plt.bar(-0.4, diff, 0.4, color=col)
        plt.ylim(0.99 * diff, 1.01*diff)
        plt.xticks([])

        sd = '{}/'.format(direct)
        # plt.savefig(sd + 'fig_{}.1.png'.format(t), dpi=256)
        # plt.close()
        plt.show()
        continue
        #
        # # PLOTS 2
        # td = '{}_finalpass/'.format(direct)
        # dst = 'rnd_{}'.format(t)
        # fs = [td + f for f in os.listdir(td) if dst in f][0]
        # cst = ChromStruct('chr1', init=fs)
        # prm, tprm = cst.params, cst.true_params
        # tclh, clh = cst.true_clh[0], cst.stat.best_lh
        #
        # plt.figure(figsize=(10, 5))
        # plt.subplots_adjust(wspace=0.055)
        # plt.subplot(1, 8, (1, 6))
        # plt.title('random params {}'.format(t))
        # x = np.arange(6)
        # ty = 10**tprm[:-1]*1e8
        # iy = 10**prm[:-1]*1e8
        # plt.bar(x-0.4, ty, 0.4, label='true params', color='black')
        # plt.bar(x, iy, 0.4, label='avg. top 3 NM', color='darkorange')
        # plt.ylabel(r'$u_{del}\ \cdot 10^{-8}$')
        # plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
        # plt.xlabel('deleterious fitness effect')
        # plt.legend()
        #
        # plt.subplot(187)
        # plt.title('pi0')
        # plt.bar(-0.4, tprm[-1], 0.4, color='black')
        # plt.bar(0, prm[-1], 0.4, color='darkorange')
        # plt.xticks([])
        # plt.yticks([])
        #
        # sp3 = plt.subplot(188)
        # sp3.yaxis.tick_right()
        # sp3.yaxis.set_label_position('right')
        # plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
        # diff = np.exp(-1e-05 * clh)-np.exp(-1e-05 * tclh)
        # col = 'green' if diff > 0 else 'red'
        # plt.bar(-0.4, diff, 0.4, color=col)
        # plt.ylim(0.99 * diff, 1.01*diff)
        # plt.xticks([])
        #
        # plt.savefig(td + 'fig_{}.2a.png'.format(t), dpi=256)
        # plt.close()

    # END


def opts(rdir):
    # opening/sorting optimization runs from octover
    fdir = root_dir + '/result/final_files/opttests/{}/'.format(rdir)

    dr = defaultdict(list)
    for f in os.listdir(fdir):
        if f.endswith('.txt'):
            fname = fdir + f
            dist = re.search('rnd_\d', fname).group()
            dr[dist].append(ChromStruct(chrom='chr1', init=fname))

    for dist in dr:
        # get the list of all runs for that distribution
        rlst = dr[dist]
        # get an indexed list of all likelihoods for run list
        ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
        # get indices of the top 3 best LH runs (after sorting on LH)
        top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

        # # determine whether the best LH is better than "true"
        # if rlst[top3[0]].stat.best_lh < rlst[top3[0]].stat.true_clh:
        #     msg = 'LH_BEST > LH_TRUE'
        # else:
        #     msg = 'LH_BEST < LH_TRUE'

        # create the plot
        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(left=0.08, wspace=0.65)
        plt.subplot(1, 9, (1, 6))
        plt.title('classic sweeps')

        x = np.arange(6)
        n = 2
        w = 0.8 / n
        s = -w / 2.0

        # get the true params for the run and plot as the first variable
        tprm = 10 ** rlst[0].stat.true_params[:-1]
        # tpmf = tprm / np.sum(tprm)
        tpmf = tprm
        plt.bar(x + s, tpmf, w, color='k', label='true')
        s += w

        # plot results from the top3 indices of runs (averaged)
        pmfs = []
        for i in top3:
            r = rlst[i]
            pmfs.append(r.avec[0])
        pmfs = np.array(pmfs)
        mpmf = np.average(pmfs, axis=0)
        plt.bar(x + s, mpmf, w, label='best')
        plt.ylabel(r'$\alpha$')
        plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
        plt.xlabel('adaptive fitness effect')
        plt.legend()

        # plot alpha/udel
        sp3 = plt.subplot(197)
        sp3.yaxis.tick_right()
        sp3.yaxis.set_label_position('right')
        plt.title(r'$\alpha$')
        plt.bar(0, sum(mpmf), 0.8)
        plt.xticks([])
        # plt.yticks([])

        # plot pi0
        plt.subplot(198)
        plt.title(r'$\pi_0$')
        pi0 = rlst[top3[0]].params[-1]
        tpi0 = rlst[top3[0]].stat.true_params[-1]
        plt.bar(0, tpi0, 1, color='k')
        plt.bar(1, pi0, 1)
        plt.xticks([])
        plt.yticks([])

        # plot clh_best - clh_true
        sp4 = plt.subplot(199)
        sp4.yaxis.tick_right()
        sp4.yaxis.set_label_position('right')
        plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
        clh, tclh = rlst[top3[0]].stat.best_lh, rlst[top3[0]].stat.true_clh
        diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
        col = 'green' if diff > 0 else 'red'
        plt.bar(-0.4, diff, 0.4, color=col)
        plt.ylim(0.99 * diff, 1.01 * diff)
        plt.xticks([])

        # create a figure for the results of the top 3 runs
        # ttl = 'INITIAL={}\nDIST={}\n{}'.format(dist, tkn, msg)
        # plt.title(ttl)
        fsave = fdir + dist + '.png'
        plt.savefig(fsave, dpi=256)
        plt.close()
        # plt.show()
        # break


def opts_2():
    # opening/sorting optimization runs from octover
    # fdir = root_dir + '/result/final_files/opttests/round_4_finalpass/'
    fdir = root_dir + '/result/final_files/opttests/round_3/test_1/'
    files = [fdir + f for f in os.listdir(fdir) if f.endswith('.txt')]

    for f in files:
        # create chromstruct from final saved file
        cst = ChromStruct(chrom='chr1', init=f)
        # use regex to get random distribution label
        dist = re.search('iprm_\d{2}', f).group()

        # create the plot
        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(left=0.08, wspace=0.65)
        plt.subplot(1, 9, (1, 6))
        plt.title('BS selection')

        x = np.arange(6)
        n = 2
        w = 0.8 / n
        s = -w / 2.0

        # get the true params for the run and plot as the first variable
        tprm = 10 ** cst.stat.true_params[:-1]
        # tpmf = tprm / np.sum(tprm)
        tpmf = tprm * 1e8
        plt.bar(x + s, tpmf, w, color='k', label='true')
        s += w

        # plot results from the top3 indices of runs (averaged)
        # pmfs = []
        # for i in top3:
        #     r = rlst[i]
        #     pmfs.append(r.avec[0])
        # pmfs = np.array(pmfs)
        # mpmf = np.average(pmfs, axis=0)
        mpmf = cst.uvec[0] * 1e8
        plt.bar(x + s, mpmf, w, label='best')
        plt.ylabel(r'$\mu_{del} \cdot 10^8$')
        plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
        plt.xlabel('deleterious fitness effect')
        plt.legend()

        # plot alpha/udel
        sp3 = plt.subplot(197)
        sp3.yaxis.tick_right()
        sp3.yaxis.set_label_position('right')
        plt.title(r'$\mu_{del} \cdot 10^8$')
        plt.bar(0, sum(mpmf), 0.8)
        plt.xticks([])
        # plt.yticks([])

        # plot pi0
        plt.subplot(198)
        plt.title(r'$\pi_0$')
        pi0 = cst.params[-1]
        tpi0 = cst.stat.true_params[-1]
        plt.bar(0, tpi0, 1, color='k')
        plt.bar(1, pi0, 1)
        plt.xticks([])
        plt.yticks([])

        # plot clh_best - clh_true
        sp4 = plt.subplot(199)
        sp4.yaxis.tick_right()
        sp4.yaxis.set_label_position('right')
        plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
        clh, tclh = cst.stat.best_lh, cst.stat.true_clh[0]
        diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
        col = 'green' if diff > 0 else 'red'
        plt.bar(-0.4, diff, 0.4, color=col)
        plt.ylim(0.99 * diff, 1.01 * diff)
        plt.xticks([])

        # create a figure for the results of the top 3 runs
        # ttl = 'INITIAL={}\nDIST={}\n{}'.format(dist, tkn, msg)
        # plt.title(ttl)
        fsave = fdir + dist + '.png'
        plt.savefig(fsave, dpi=256)
        plt.close()
        # plt.show()
        # break


def opts_3():
    # opening/sorting optimization runs from octover
    fdir = root_dir + '/result/final_files/opttests/bs1cs1/'

    dr = defaultdict(list)
    for f in os.listdir(fdir):
        if f.endswith('.txt'):
            fname = fdir + f
            dist = re.search('rnd_\d', fname).group()
            dr[dist].append(ChromStruct(chrom='chr1', init=fname))

    for dist in dr:
        # get the list of all runs for that distribution
        rlst = dr[dist]
        # get an indexed list of all likelihoods for run list
        ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
        # get indices of the top 3 best LH runs (after sorting on LH)
        top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

        # # determine whether the best LH is better than "true"
        # if rlst[top3[0]].stat.best_lh < rlst[top3[0]].stat.true_clh:
        #     msg = 'LH_BEST > LH_TRUE'
        # else:
        #     msg = 'LH_BEST < LH_TRUE'

        # create the plot
        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(left=0.08, wspace=0.65, hspace=0.65)

        # panel 1: BS
        plt.subplot(2, 9, (1, 6))
        plt.title('background selection')

        x = np.arange(6)
        n = 2
        w = 0.8 / n
        s = -w / 2.0

        # get the true params for the run and plot as the first variable
        tprm = 10 ** rlst[0].stat.true_params[0:6]
        tpmf = tprm * 1e8
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

        # plot udel;
        sp3 = plt.subplot(2, 9, 7)
        sp3.yaxis.tick_right()
        sp3.yaxis.set_label_position('right')
        plt.title(r'$\mu_{del} \cdot 10^8$')
        plt.bar(0, sum(mpmf), 0.8)
        plt.xticks([])


        # panel 2: CS
        plt.subplot(2, 9, (10, 15))
        plt.title('classic sweeps')
        x = np.arange(6)
        n = 2
        w = 0.8 / n
        s = -w / 2.0

        # get the true params for the run and plot as the first variable
        tprm = 10 ** rlst[0].stat.true_params[6:12]
        tpmf = tprm
        plt.bar(x + s, tpmf, w, color='k', label='true')
        s += w

        # plot results from the top3 indices of runs (averaged)
        pmfs = []
        for i in top3:
            r = rlst[i]
            pmfs.append(r.avec[0])
        pmfs = np.array(pmfs)
        mpmf = np.average(pmfs, axis=0)
        plt.bar(x + s, mpmf, w, label='best')
        plt.ylabel(r'$\alpha$')
        plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
        plt.xlabel('adaptive fitness effect')
        plt.legend()

        # plot alpha
        sp3 = plt.subplot(2, 9, 16)
        sp3.yaxis.tick_right()
        sp3.yaxis.set_label_position('right')
        plt.title(r'$\alpha$')
        plt.bar(0, sum(mpmf), 0.8)
        plt.xticks([])

        # plot pi0
        plt.subplot(198)
        plt.title(r'$\pi_0$')
        pi0 = rlst[top3[0]].params[-1]
        tpi0 = rlst[top3[0]].stat.true_params[-1]
        plt.bar(0, tpi0, 1, color='k')
        plt.bar(1, pi0, 1)
        plt.xticks([])
        plt.yticks([])

        # plot clh_best - clh_true
        sp4 = plt.subplot(199)
        sp4.yaxis.tick_right()
        sp4.yaxis.set_label_position('right')
        plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
        clh, tclh = rlst[top3[0]].stat.best_lh, rlst[top3[0]].stat.true_clh
        diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
        col = 'green' if diff > 0 else 'red'
        plt.bar(-0.4, diff, 0.4, color=col)
        plt.ylim(0.99 * diff, 1.01 * diff)
        plt.xticks([])

        # create a figure for the results of the top 3 runs
        # ttl = 'INITIAL={}\nDIST={}\n{}'.format(dist, tkn, msg)
        # plt.title(ttl)
        fsave = fdir + dist + '.png'
        plt.savefig(fsave, dpi=256)
        plt.close()
        # plt.show()
        # break


def opts_4(rdir):
    # opening/sorting optimization runs from octover
    fdir = root_dir + '/result/final_files/bound_bsx_feb2019/'

    dr = defaultdict(list)
    for f in os.listdir(fdir):
        if f.endswith('.txt'):
            fname = fdir + f
            dist = re.search('rnd_\d', fname).group()
            dr[dist].append(ChromStruct(chrom='chr1', init=fname))

    for dist in dr:
        # get the list of all runs for that distribution
        rlst = dr[dist]
        # get an indexed list of all likelihoods for run list
        ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
        # get indices of the top 3 best LH runs (after sorting on LH)
        top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

        # create the plot
        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(left=0.08, wspace=0.65)
        plt.subplot(1, 9, (1, 6))
        plt.title('background selection')

        x = np.arange(6)
        n = 2
        w = 0.8 / n
        s = -w / 2.0

        # get the true params for the run and plot as the first variable
        tprm = 10 ** rlst[0].stat.true_params[:-1]
        tpmf = tprm * 1e8
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
        sp3 = plt.subplot(197)
        sp3.yaxis.tick_right()
        sp3.yaxis.set_label_position('right')
        plt.title(r'$\mu_{del} \cdot 10^8$')
        plt.bar(0, sum(mpmf), 0.8)
        plt.xticks([])

        # plot pi0
        plt.subplot(198)
        plt.title(r'$\pi_0$')
        pi0 = rlst[top3[0]].params[-1]
        tpi0 = rlst[top3[0]].stat.true_params[-1]
        plt.bar(0, tpi0, 1, color='k')
        plt.bar(1, pi0, 1)
        plt.xticks([])
        plt.yticks([])

        # plot clh_best - clh_true
        sp4 = plt.subplot(199)
        sp4.yaxis.tick_right()
        sp4.yaxis.set_label_position('right')
        plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
        clh, tclh = rlst[top3[0]].stat.best_lh, rlst[top3[0]].stat.true_clh
        diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
        col = 'green' if diff > 0 else 'red'
        plt.bar(-0.4, diff, 0.4, color=col)
        plt.ylim(0.99 * diff, 1.01 * diff)
        plt.xticks([])

        # create a figure for the results of the top 3 runs
        # ttl = 'INITIAL={}\nDIST={}\n{}'.format(dist, tkn, msg)
        # plt.title(ttl)
        fsave = fdir + dist + '.png'
        plt.savefig(fsave, dpi=256)
        plt.close()
        # plt.show()
        # break

# opts('cs1_smpl_2')
opts_4('bs1_detrm')
# some_plots()