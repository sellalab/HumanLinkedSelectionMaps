__author__ = 'davidmurphy'


import os
import seaborn
import tarfile
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir, \
    human_autosomes, RunStruct, cst_from_fldr
from classes.phylotree import parse_exptotsub
from data_processing.functions import rsquared_function

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


#%% BASIC SORT (MORE BINS)


def basic_sort_plot(fldr, num, label, subtr=False):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    sdir = root_dir + '/result/final_files/sfigs/'
    sort_file = fdir + 'basic_sort_n{}.txt'.format(num)
    div, pi, pred = np.loadtxt(sort_file).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    # normalize by pi0
    pi /= pi0
    pred /= pi0

    # create new plot
    plt.figure(figsize=(2.5, 2.5))
    plt.subplots_adjust(top=0.98, right=0.95, left=0.2, bottom=0.17)

    # subtract pred-pi plot
    if subtr:
        plt.plot(pred, pi-pred, marker='o', ms=5, markerfacecolor='None',
                 markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
                 alpha=0.75)
        # ytick = np.arange(0.5, 1.3, 0.1)
        # plt.yticks(ytick)
        plt.ylabel(r'observed - predicted $\pi/\pi_0$', labelpad=3)

        plt.plot([0.5, 1.02], [0, 0], color='darkslategray',
                 ls='--', alpha=0.65)

        plt.text(0.76, 0.3, r'$n_{bins}=$' + str(num), ha='center',
                 va='center')
        plt.yticks(np.arange(-0.15, 0.4, 0.05))
        plt.ylim(-0.18, 0.38)
        plt.legend(loc='lower right', ncol=3)
        f_save = sdir + '{}_n{}_sortplot_subtr.png'.format(fldr, num)

    # standard plot
    else:
        rsq = rsquared_function(pi, pred)
        print rsq
        plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
                 markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
                 alpha=0.75, label=label)
        ytick = np.arange(0.45, 1.3, 0.1)
        plt.yticks(ytick)
        plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)

        plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
                 ls='--', alpha=0.65)
        # plt.text(0.52, 1.17, r'$n_{bins}=$' + str(num))
        plt.axhline(y=1, color='k', alpha=0.8, ls='-')
        plt.text(0.72, 1.02, 'without linked selection', ha='center',
                 va='center')
        plt.ylim(0.45, 1.15)
        plt.legend(loc='lower right', ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', prop=dict(size=7))
        f_save = sdir + '{}_n{}_sortplot.png'.format(fldr, num)
        plt.text(0.25, 0.6, r'$R^2$={}'.format(round(rsq, 4)),
                 transform=plt.gcf().transFigure)

    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    xtick = np.arange(0.25, 1.2, 0.1)
    plt.xticks(xtick)
    plt.xlim(0.45, 1.15)

    plt.savefig(f_save, dpi=512)


# for n in [100, 250, 500, 1000, 2000]:
# for i in range(1, 21):
mnb_tkn = '603 520 444 375'.split()
bth_tkn = '500 550 600 650'.split()
# lab = r'precalc $B \geq 0.2$'
lab = r'LLH $B \geq 0.2$'
for i in range(len(mnb_tkn)):
    lab = r'CADD LLH $B \geq{}$'.format(0.001*float(bth_tkn[i]))
    fldr = 'cadd93_bth{}'.format(bth_tkn[i])
    basic_sort_plot(fldr, 100, lab, subtr=False)
for i in range(len(mnb_tkn)):
    lab = r'CADD prec $B \geq{}$'.format(0.001*float(bth_tkn[i]))
    fldr = 'cadd93_mnb{}'.format(mnb_tkn[i])
    basic_sort_plot(fldr, 100, lab, subtr=False)
for i in range(len(mnb_tkn)):
    lab = r'ape LLH $B \geq{}$'.format(0.001*float(bth_tkn[i]))
    fldr = 'ape_cons94_bth{}'.format(bth_tkn[i])
    basic_sort_plot(fldr, 100, lab, subtr=False)
for i in range(len(mnb_tkn)):
    lab = r'ape prec $B \geq{}$'.format(0.001*float(bth_tkn[i]))
    fldr = 'ape_cons94_mnb{}'.format(mnb_tkn[i])
    basic_sort_plot(fldr, 100, lab, subtr=False)
#%% SORT PLOT DETAIL CIRCA B=1

def sort_plot_detail(fldr, num, pct=0.98):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    sdir = root_dir + '/result/final_files/sfigs/'
    sort_file = fdir + 'basic_sort_n{}.txt'.format(num)
    div, pi, pred = np.loadtxt(sort_file).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    # normalize by pi0
    pi /= pi0
    pred /= pi0

    # get the index corresponding to the percentile of sites to take
    idx = int(pct * num)
    # take only the top percentile of sites
    pi = pi[idx:]
    pred = pred[idx:]

    # create new plot
    plt.figure(figsize=(2.5, 2.5))
    plt.subplots_adjust(top=1, right=1, left=0.2, bottom=0.17)

    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
             alpha=0.75)
    # ytick = np.arange(0.5, 1.3, 0.1)
    # plt.yticks(ytick)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)

    # plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
    #          ls='--', alpha=0.65)
    # plt.text(0.52, 1.17, r'$n_{bins}=$' + str(num))
    plt.axhline(y=1, color='k', alpha=0.8, ls='-',
                label='without linked selection')
    # plt.text(0.775, 1.02, 'without linked selection', ha='center',
    #          va='center')
    # plt.ylim(0.5, 1.12)
    plt.legend(loc='upper left', ncol=3)
    f_save = sdir + '{}_n{}_sortplot_detail.png'.format(fldr, num)

    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    # xtick = np.arange(0.5, 1.01, 0.1)
    # plt.xticks(xtick)
    # plt.xlim(0.5, 1.02)

    plt.savefig(f_save, dpi=512)


sort_plot_detail('cadd93', 2000)
#%% ADDITIONAL SORTED PLOTS


def sort_more_data(fldr, titl):
    # load results
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    if 'standard sorting' in titl:
        f_sort = fdir + 'predsort_more_stdsort.txt'
        f_name = '/{}_sort_more_stdsort.png'.format(fldr)
    else:
        f_sort = fdir + 'predsort_more.txt'
        f_name = '/{}_sort_more.png'.format(fldr)

    div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    # normalize by pi0
    pi /= pi0
    pred /= pi0

    # create new plot
    plt.figure(figsize=(4, 12))
    plt.subplots_adjust(top=0.97, right=0.99, left=0.2, bottom=0.04)
    # if 'clean' in fldr:
    #     spec = fldr.split('_')[0]
    #     spec = spec[0].upper() + spec[1:]
    #     plt.suptitle('{} 6% Most Conserved'.format(spec))
    # else:
    #     plt.suptitle('CADD Top 7%')
    plt.suptitle(titl, y=0.99)
    # plot standard predicted/observed on top
    plt.subplot(411)
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
             alpha=0.75, label=r'$\pi/\pi_0$')
    plt.plot(pred, div, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='crimson', markeredgewidth=0.9, lw=0,
             alpha=0.75, label=r'$D/\bar{D}$')
    xtick = np.arange(0.5, 1.01, 0.1)
    ytick = np.arange(0.4, 1.3, 0.2)
    plt.xticks(xtick, color='white')
    # plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    # plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.ylabel('observed values', labelpad=3)

    plt.yticks(ytick)
    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.ylim(0.35, 1.35)
    plt.xlim(0.55, 1.02)
    plt.legend(loc='lower right', ncol=3)

    # plot cM/Mb
    plt.subplot(412)
    plt.plot(pred, np.log10(cm * 1e6), marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='deepskyblue', markeredgewidth=0.9, lw=0,
             alpha=0.75)
    plt.xticks(xtick, color='white')
    plt.xlim(0.55, 1.02)
    plt.ylabel('recombination rate\n(log10 cM/Mb)', labelpad=3)
    plt.ylim(-1.2, 1.2)

    # plot conservation
    plt.subplot(413)
    plt.plot(pred, cn * 1e-3, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkviolet', markeredgewidth=0.9, lw=0,
             alpha=0.75)
    plt.xticks(xtick, color='white')
    # plt.yticks(np.arange(0.38, 0.44, 0.01))
    plt.xlim(0.55, 1.02)
    plt.ylim(0.01, 0.13)
    plt.ylabel('conservation\nscore', labelpad=3)

    # plot GC content
    plt.subplot(414)
    plt.plot(pred, gc, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='forestgreen', markeredgewidth=0.9, lw=0,
             alpha=0.75)
    plt.xticks(xtick, y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.yticks(np.arange(0.38, 0.47, 0.02))
    plt.xlim(0.55, 1.02)
    plt.ylim(0.37, 0.462)
    plt.ylabel('GC fraction', labelpad=3)

    sdir = root_dir + '/result/final_files/sfigs'
    f_save = sdir + f_name

    plt.savefig(f_save, dpi=256)
    plt.close()


#%% SORT DIVERGENCE DATA


def sort_divergence(fldr):
    """sort additional divergence data for B=1 analyses"""
    # load results
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'predsort_more.txt'
    f_bbin = fdir + 'bbin_rates.txt'
    f_save = '/{}_bbin_sort.png'.format(fldr)

    bbin, at_rate, a2t_rate, gc_rate, all_rate, hm_rate = np.loadtxt(f_bbin).T
    div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    norm_all_rate = all_rate / np.mean(all_rate)
    norm_at_rate = at_rate / np.mean(at_rate)
    norm_gc_rate = gc_rate / np.mean(gc_rate)
    norm_a2t_rate = a2t_rate / np.mean(a2t_rate)
    norm_hm_rate = hm_rate / np.mean(hm_rate)

    # normalize by pi0
    pi /= pi0
    pi *= div
    pi /= norm_all_rate
    pred /= pi0

    # create new plot
    plt.figure(figsize=(3.5, 6))
    plt.subplots_adjust(top=0.99, right=0.98, left=0.16, bottom=0.08,
                        hspace=0.08)

    # plot standard predicted/observed on top
    plt.subplot(211)
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
             alpha=0.75, label='CADD 7%')
    xtick = np.arange(0.6, 1.01, 0.1)
    ytick = np.arange(0.6, 1.11, 0.1)
    plt.xticks(xtick, color='white')
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick)
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.ylim(0.55, 1.15)
    plt.xlim(0.55, 1.02)
    plt.legend(loc='upper left', ncol=1)

    # plot divergence rates
    plt.subplot(212)
    plt.plot(pred, norm_all_rate, marker='o', ms=5, color='black', lw=0,
             alpha=0.75, label='all')
    plt.plot(pred, norm_at_rate, marker='o', ms=5, color='dodgerblue', lw=0,
             alpha=0.75, label='AT')
    plt.plot(pred, norm_gc_rate, marker='o', ms=5, color='fuchsia', lw=0,
             alpha=0.75, label='GC')
    # plt.plot(pred, norm_a2t_rate, marker='o', ms=5, color='purple', lw=0,
    #          alpha=0.75, label='A>T')
    plt.plot(pred, norm_hm_rate, marker='o', ms=5, color='darkturquoise', lw=0,
             alpha=0.75, label='HM')
    plt.xlim(0.55, 1.02)
    plt.legend(loc='upper left', ncol=2)

    sdir = root_dir + '/result/final_files/sfigs'
    f_save = sdir + f_save
    plt.xticks(xtick)
    plt.ylabel(r'$D/\bar{D}$', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.savefig(f_save, dpi=512)
    plt.close()


sort_divergence('cadd93')


#%% SORT DIVERGENCE FOR BINNED VS. PHYLOFIT VERSIONS


def sort_divergence_2(fldr):
    """sort additional divergence data for B=1 analyses"""
    # load results
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'predsort_more.txt'
    f_bbin = fdir.replace('_new_nu', '') + 'bbin_rates.txt'
    f_save = '/{}_bbin_sort_2.png'.format(fldr)

    bbin, at_rate, a2t_rate, gc_rate, all_rate, hm_rate = np.loadtxt(f_bbin).T
    div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    norm_all_rate = all_rate / np.mean(all_rate)

    # normalize by pi0
    pi /= pi0
    pi *= div
    pi /= norm_all_rate
    pred /= pi0

    # create new plot
    plt.figure(figsize=(3.5, 6))
    plt.subplots_adjust(top=0.99, right=0.98, left=0.16, bottom=0.08,
                        hspace=0.08)

    # plot standard predicted/observed on top
    plt.subplot(211)
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
             alpha=0.75, label='CADD 7%')
    xtick = np.arange(0.6, 1.01, 0.1)
    ytick = np.arange(0.6, 1.11, 0.1)
    plt.xticks(xtick, color='white')
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick)
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.ylim(0.55, 1.15)
    plt.xlim(0.55, 1.02)
    plt.legend(loc='upper left', ncol=1)

    # plot divergence rates
    plt.subplot(212)
    plt.plot(pred, norm_all_rate, marker='o', ms=5, color='black', lw=0,
             alpha=0.75, label='phyloFit D')
    plt.plot(pred, div, marker='o', ms=5, color='dodgerblue', lw=0,
             alpha=0.75, label='binned D')

    plt.xlim(0.55, 1.02)
    plt.legend(loc='upper left', ncol=2)

    sdir = root_dir + '/result/final_files/sfigs'
    f_save = sdir + f_save
    plt.xticks(xtick)
    plt.ylabel(r'$D/\bar{D}$', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.savefig(f_save, dpi=512)
    plt.close()


sort_divergence_2('cadd93')


#%% SORT DIVERGENCE FOR EACH MUTATION TYPE


def sort_divergence_3(fldr):
    """sort additional divergence data for B=1 analyses"""
    # load results
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'predsort_more.txt'
    f_bbin = fdir.replace('_new_nu', '') + 'bbin_rates.txt'
    f_save = '/{}_bbin_sort_3.png'.format(fldr)

    # set directories used in run
    bbin_dir = root_dir + '/data/phast/bbins'
    res = []
    for bbin in range(100):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        f_name = '{}/bbin{}.exptotsub'.format(bbin_dir, bbin)
        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m
        # adjust total for number of branches
        for i in range(4):
            mat[i,i] /= 14.0
        # get A>C/T>G
        at_tot = (mat[0, :].sum() + mat[3, :].sum())
        actg = (mat[0, 1] + mat[3, 2]) / at_tot
        # get A>G/T>C
        agtc = (mat[0, 2] + mat[3, 1]) / at_tot
        # get A>T/T>A
        atta = (mat[0, 3] + mat[3, 0]) / at_tot
        # get C>A/G>T
        cg_tot = (mat[1, :].sum() + mat[2, :].sum())
        cagt = (mat[1, 0] + mat[2, 3]) / cg_tot
        # get C>G/G>C
        cggc = (mat[1, 2] + mat[2, 1]) / cg_tot
        # get C>T/G>A
        ctga = (mat[1, 3] + mat[2, 0]) / cg_tot

        row = actg, agtc, atta, cagt, cggc, ctga
        res.append(row)

    res = np.array(res)
    # actg, agtc, atta, cagt, cggc, ctga = res.T
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown', 'darkturquoise']

    bbin, at_rate, a2t_rate, gc_rate, all_rate, hm_rate = np.loadtxt(f_bbin).T
    div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    norm_all_rate = all_rate / np.mean(all_rate)

    # normalize by pi0
    pi /= pi0
    pi *= div
    pi /= norm_all_rate
    pred /= pi0

    # create new plot
    plt.figure(figsize=(3.5, 6))
    plt.subplots_adjust(top=0.99, right=0.98, left=0.17, bottom=0.08,
                        hspace=0.08)

    # plot standard predicted/observed on top
    plt.subplot(211)
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
             alpha=0.75, label='CADD 7%')
    xtick = np.arange(0.6, 1.01, 0.1)
    ytick = np.arange(0.6, 1.11, 0.1)
    plt.xticks(xtick, color='white')
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick)
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.ylim(0.55, 1.15)
    plt.xlim(0.55, 1.02)
    plt.legend(loc='upper left', ncol=1)

    # plot divergence rates
    plt.subplot(212)
    for (i, r) in enumerate(res.T):
        plt.plot(pred, r/r.mean(), marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])

    plt.xlim(0.55, 1.02)
    plt.legend(loc='upper left', ncol=2)

    sdir = root_dir + '/result/final_files/sfigs'
    f_save = sdir + f_save
    plt.xticks(xtick)
    plt.ylabel(r'$D/\bar{D}$', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.savefig(f_save, dpi=512)
    plt.close()


sort_divergence_3('cadd93')


#%% SORT DATA BY B VALUES (PI, CMMB, CONS, GC)


def sort_data(fldr, nbins=100):
    """sort additional divergence data for B=1 analyses"""
    # load results
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_n{}.txt'.format(nbins)
    # f_bbin = fdir + 'bbin_rates.txt'
    f_save1 = '/{}_sort_pi_cmmb_n{}.png'.format(fldr, nbins)
    f_save2 = '/{}_sort_cons_gc_n{}.png'.format(fldr, nbins)

    # bbin, at_rate, a2t_rate, gc_rate, all_rate, hm_rate = np.loadtxt(f_bbin).T
    div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    # div2, pi2, pred2 = np.loadtxt(fdir + '/basic_sort_n100.txt').T

    if fldr == 'cadd93_extel_rm_CG':
        rst = cst_from_fldr('cadd93_bth600')
    else:
        rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    # pimean = pi.mean()
    # norm_all_rate = all_rate / np.mean(all_rate)

    # normalize by pi0
    # pi *= (rst.stat.meanpi / (pimean*pi0))
    # fix divergence correction with best result
    # pi *= div
    # pi /= norm_all_rate
    # normalize prediction by pi0
    pi /= pi0
    pred /= pi0
    # pi2 /= pi0
    # pred2 /= pi0

    # create new plot
    plt.figure(figsize=(3, 21.0/4.0))
    plt.subplots_adjust(top=0.99, right=0.98, left=0.22, bottom=0.08,
                        hspace=0.08)

    # plot standard predicted/observed on top
    plt.subplot(211)
    plt.plot(pred, pi, marker='o', ms=5, lw=0, alpha=0.75,
             color='darkorange', label='CADD 7%')
    # plt.plot(pred2, pi2, marker='X', ms=5, lw=0, alpha=0.75,
    #          color='k')
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    xtick = np.arange(0.5, 1.01, 0.1)
    ytick = np.arange(0.5, 1.41, 0.1)
    plt.xticks(xtick, color='white')
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.yticks(ytick)
    plt.ylim(0.52, 1.15)
    plt.xlim(0.52, 1.02)
    plt.legend(loc='lower right')
    plt.text(0.23, 0.96, 'A', transform=plt.gcf().transFigure)

    # plot divergence rates
    plt.subplot(212)
    plt.plot(pred, np.log10(cm * 1e6), marker='o', ms=5, lw=0, alpha=0.75,
             color='deepskyblue', label='CADD 7%')

    plt.xticks(xtick)
    plt.ylabel('recombination rate (log10 cM/Mb)', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylim(-1.2, 1.2)
    plt.xlim(0.52, 1.02)
    plt.text(0.23, 0.485, 'B', transform=plt.gcf().transFigure)

    sdir = root_dir + '/result/final_files/sfigs'
    f_save = sdir + f_save1
    plt.savefig(f_save, dpi=512)
    plt.close()

    # create new plot
    plt.figure(figsize=(3, 21.0/4.0))
    plt.subplots_adjust(top=0.99, right=0.98, left=0.22, bottom=0.08,
                        hspace=0.08)

    # plot standard predicted/observed on top
    plt.subplot(211)
    plt.plot(pred, cn/1000.0, marker='o', ms=5, lw=0, alpha=0.75,
             color='purple')
    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick, color='white')
    plt.ylabel('conservation score', labelpad=3)
    plt.xlim(0.52, 1.02)
    plt.legend(loc='upper left')
    plt.text(0.23, 0.96, 'C', transform=plt.gcf().transFigure)

    # plot divergence rates
    plt.subplot(212)
    plt.plot(pred, gc, marker='o', ms=5, lw=0, alpha=0.75,
             color='forestgreen')
    plt.xticks(xtick)
    plt.ylabel('GC fraction', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.xlim(0.52, 1.02)
    plt.text(0.23, 0.485, 'D', transform=plt.gcf().transFigure)

    sdir = root_dir + '/result/final_files/sfigs'
    f_save = sdir + f_save2
    plt.savefig(f_save, dpi=512)
    plt.close()


sort_data('cadd93_extel_rm_CG', nbins=100)


#%% LOAD SUBSTITUTION RATES AND SNP TYPES


def get_sub_rates(fldr):
    rst = cst_from_fldr(fldr)
    bbin_dir = root_dir + '/data/phast/bbins'
    res = []
    for bbin in range(100):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        # f_name = '{}/bbin{}.exptotsub'.format(bbin_dir, bbin)
        f_name = '{}/{}.bbin{}.exptotsub'.format(bbin_dir, rst.tkn, bbin)

        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m
        # adjust total for number of branches
        for i in range(4):
            mat[i,i] /= 14.0
        # get A>C/T>G
        at_tot = (mat[0, :].sum() + mat[3, :].sum())
        # at_tot = (mat[0, 0].sum() + mat[3, 3].sum())

        actg = (mat[0, 1] + mat[3, 2]) / at_tot
        # get A>G/T>C
        agtc = (mat[0, 2] + mat[3, 1]) / at_tot
        # get A>T/T>A
        atta = (mat[0, 3] + mat[3, 0]) / at_tot
        # get C>A/G>T
        cg_tot = (mat[1, :].sum() + mat[2, :].sum())
        # cg_tot = (mat[1, 1].sum() + mat[2, 2].sum())
        cagt = (mat[1, 0] + mat[2, 3]) / cg_tot
        # get C>G/G>C
        cggc = (mat[1, 2] + mat[2, 1]) / cg_tot
        # get C>T/G>A
        ctga = (mat[1, 3] + mat[2, 0]) / cg_tot

        row = actg, agtc, atta, cagt, cggc, ctga, at_tot, cg_tot
        res.append(row)

    res = np.array(res)

    dv_vals = []
    for r in res[:,:6].T:
        dv = r/r.mean()
        dv_vals.append(dv)
    dv = np.array(dv_vals)

    return res, dv


# fldr = 'cadd93_bth600'
# res, dv = get_sub_rates(fldr)


#%%
def get_sub_rates_2(fldr):
    rst = cst_from_fldr(fldr)
    bbin_dir = root_dir + '/data/phast/bbins'
    res = []
    msk = []
    for bbin in range(1500, 2000):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        # f_name = '{}/bbin{}.exptotsub'.format(bbin_dir, bbin)
        f_name = '{}/{}.bbin{}.exptotsub'.format(bbin_dir, rst.tkn, bbin)
        if not os.path.isfile(f_name):
            print '{} missing bbin {}'.format(fldr, bbin)
            msk.append(False)
            continue

        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m
        # adjust total for number of branches
        for i in range(4):
            mat[i,i] /= 14.0
        # get A>C/T>G
        at_tot = (mat[0, :].sum() + mat[3, :].sum())
        # at_tot = (mat[0, 0].sum() + mat[3, 3].sum())

        actg = (mat[0, 1] + mat[3, 2]) / at_tot
        # get A>G/T>C
        agtc = (mat[0, 2] + mat[3, 1]) / at_tot
        # get A>T/T>A
        atta = (mat[0, 3] + mat[3, 0]) / at_tot
        # get C>A/G>T
        cg_tot = (mat[1, :].sum() + mat[2, :].sum())
        # cg_tot = (mat[1, 1].sum() + mat[2, 2].sum())
        cagt = (mat[1, 0] + mat[2, 3]) / cg_tot
        # get C>G/G>C
        cggc = (mat[1, 2] + mat[2, 1]) / cg_tot
        # get C>T/G>A
        ctga = (mat[1, 3] + mat[2, 0]) / cg_tot

        row = actg, agtc, atta, cagt, cggc, ctga, at_tot, cg_tot
        res.append(row)
        msk.append(True)

    res = np.array(res)

    dv_vals = []
    for r in res[:,:6].T:
        dv = r/r.mean()
        dv_vals.append(dv)
    dv = np.array(dv_vals)

    return res, dv, np.array(msk)


fldr = 'cadd93_extel_rm_CG'
res, dv, msk = get_sub_rates_2(fldr)


# LOAD *FIXED* SNP DATA


def get_snp_poly(fldr, nbins=100):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    # f_sort = fdir + 'predsort_more.txt'
    # f_sort = fdir + 'basic_sort_n{}.txt'.format(nbins)
    f_sort = fdir + 'sort_gc_cm_cn_n{}.txt'.format(nbins)

    f_ssnp = fdir + 'sort_snptype_fix_anc_n{}.txt'.format(nbins)

    snp = np.loadtxt(f_ssnp)
    # div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    div, pi, pred = np.loadtxt(f_sort)[:,:3].T

    if fldr == 'cadd93_extel_rm_CG':
        rst = cst_from_fldr('cadd93_bth600')
    else:
        rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    all_pimean = pi.mean()
    # pred_mean = pred.mean()

    # pred *= (all_pimean/(pred_mean*pi0))
    pred /= pi0
    pi_vals = []
    for s in snp[:,:6].T:
        pi = s*all_pimean/(s.mean()*pi0)
        pi_vals.append(pi)
    ancpi = np.array(pi_vals)

    return ancpi, pred, snp


# fldr = 'cadd93_bth600'
ancpi, pred, snp = get_snp_poly(fldr, nbins=2000)
#%%
ancpi = ancpi.T[1500:][msk].T
pred = pred[1500:][msk]
snp = snp[1500:][msk]

#%% AT ORIGINATING 2
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
f_save = root_dir + '/result/final_files/sfigs/{}_AT_2000_rm_CG.png'.format(fldr)
plt.figure(figsize=(3, 21.0 / 4.0))
plt.subplots_adjust(top=0.95, right=0.975, hspace=0.1, bottom=0.1, left=0.18)
plt.subplot(211)
plt.title('AT originating mutations')
for i in [0, 1, 2]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
# xtick = np.arange(0.5, 1.01, 0.1)
# ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
# plt.yticks(ytick)
plt.ylim(0.6, 1.8)
plt.xlim(0.9, 1.0)
plt.legend(loc='upper left', ncol=2, frameon=1, framealpha=0.75,
           facecolor='white')
# plt.text(0.01, 0.98, 'A', transform=plt.gcf().transFigure)

plt.subplot(212)
for i in [0, 1, 2]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
# xtick = np.arange(0.5, 1.01, 0.1)
# plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.9, 1.0)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
# ytick = np.arange(0.9, 1.31, 0.1)
# plt.yticks(ytick)
# plt.ylim(0.82, 1.38)
plt.legend(loc='upper left', frameon=1, framealpha=0.75, facecolor='white')
# plt.text(0.01, 0.47, 'B', transform=plt.gcf().transFigure)
plt.savefig(f_save, dpi=512)
plt.close()
#%% UN-SCALED SUBSTITUTION RATES
def unscaled_divergence(fldr, res, pred):
    f_save = final_dir + '/sfigs/fig_S30B.{}.png'.format(fldr)
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
    plt.figure(figsize=(3,3))
    plt.subplots_adjust(top=0.94, right=0.99, hspace=0.1, bottom=0.12, left=0.16)
    plt.title('Substitution rates')
    for (i,r) in enumerate(res[:,:6].T):
        plt.plot(pred, r, marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel('substitution rate', labelpad=3)
    plt.xlim(0.58, 1.02)
    plt.legend(loc='upper left', ncol=2)
    plt.text(0.01, 0.96, 'B', transform=plt.gcf().transFigure)
    plt.savefig(f_save, dpi=512)
    plt.close()

fldr = 'cadd93_bth600'
unscaled_divergence(fldr, res, pred)


#%% UN-SCALED HETEROZYGOSITY BY SNP TYPE (Fig. S30)
def unscaled_diversity(fldr, snp, pred):
    f_save = final_dir + '/sfigs/fig_S30A.{}.png'.format(fldr)
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    plt.figure(figsize=(3.2, 3.2))
    plt.subplots_adjust(top=0.92, right=0.99, hspace=0.1, bottom=0.15, left=0.2)
    plt.title('Heterozygosity by SNP')
    for (i,s) in enumerate(snp[:,:6].T):
        plt.plot(pred, s*1e3, marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    xtick = np.arange(0.6, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel(r'observed $\pi$ ($\times 10^3$)', labelpad=3)
    plt.xlim(0.58, 1.02)
    plt.legend(loc='upper left', ncol=2)
    plt.text(0.01, 0.96, 'A', transform=plt.gcf().transFigure)
    plt.savefig(f_save, dpi=512)
    plt.close()

fldr = 'cadd93_bth600'
unscaled_diversity(fldr, snp, pred)


#%% BGC ACCELERATED/DECELERATED
f_save = root_dir + '/result/final_files/sfigs/BGC_group_3.png'
labs = ['A>C/T>G (decelerated)', 'A>G/T>C (decelerated)', 'A>T/T>A',
        'C>A/G>T (accelerated)', 'C>G/G>C', 'C>T/G>A (accelerated)']
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']

plt.figure(figsize=(4, 7))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('BGC associated mutations')
for i in [0, 1, 3, 5]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
# plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.42)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')

plt.subplot(212)
for i in [0, 1, 3, 5]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
plt.legend(loc='upper left')
plt.savefig(f_save, dpi=512)
plt.close()

# non-BGC
f_save = root_dir + '/result/final_files/sfigs/non_BGC_group_3.png'
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
plt.figure(figsize=(4, 7))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('non-BGC associated mutations')
for i in [2, 4]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
# plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.42)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')

plt.subplot(212)
for i in [2, 4]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.legend(loc='upper left')
plt.savefig(f_save, dpi=512)
plt.close()

#%% TRANSITIONS / TRANSVERSIONS
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C',
        'C>T/G>A (includes CpG transitions)']
f_save = root_dir + '/result/final_files/sfigs/transitions_group_3.png'
plt.figure(figsize=(4, 7))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('transition mutations')
for i in [1, 5]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.42)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')
plt.subplot(212)
for i in [1, 5]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.legend(loc='upper left')
plt.savefig(f_save, dpi=512)
plt.close()

#TRANSVERSIONS
# ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
f_save = root_dir + '/result/final_files/sfigs/transversions_group_3.png'
plt.figure(figsize=(4, 7))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('transversion mutations')
for i in [0, 2, 3, 4]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.42)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')
plt.subplot(212)
for i in [0, 2, 3, 4]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.legend(loc='upper left')
plt.savefig(f_save, dpi=512)
plt.close()

#%% AT ORIGINATING
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
f_save = root_dir + '/result/final_files/sfigs/fig_S29.AT_{}.png'.format(fldr)
plt.figure(figsize=(3, 21.0 / 4.0))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('AT originating mutations')
for i in [0, 1, 3]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.27)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')
plt.text(0.01, 0.98, 'A', transform=plt.gcf().transFigure)

plt.subplot(212)
for i in [0, 1]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.legend(loc='upper left')
plt.text(0.01, 0.47, 'B', transform=plt.gcf().transFigure)
plt.savefig(f_save, dpi=512)
plt.close()

# GC ORIGINATING
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
f_save = root_dir + '/result/final_files/sfigs/fig_S29.GC_{}.png'.format(fldr)
plt.figure(figsize=(3, 21.0 / 4.0))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('GC originating mutations')
for i in [3, 4, 5]:
    plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.27)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')
plt.text(0.01, 0.98, 'C', transform=plt.gcf().transFigure)

plt.subplot(212)
for i in [3, 4, 5]:
    plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.legend(loc='upper left')
plt.text(0.01, 0.47, 'D', transform=plt.gcf().transFigure)
plt.savefig(f_save, dpi=512)
plt.close()


#%% MUTATIONS COMBINED
f_save = root_dir + '/result/final_files/sfigs/combined_groups_3.png'
plt.figure(figsize=(4, 7))
plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
plt.subplot(211)
plt.title('Combined mutations')
clist = ['darkturquoise', 'salmon', 'purple', 'steelblue']
labs = ['A>C/C>A/T>G/G>T', 'A>G/G>A/T>C/C>T', 'A>T/T>A', 'C>G/G>C']
pi1 = (ancpi[0] + ancpi[3]) / 2.0
pi2 = (ancpi[1] + ancpi[5]) / 2.0
pi3 = ancpi[2]
pi4 = ancpi[4]
dv1 = (dv[0] + dv[3]) / 2.0
dv2 = (dv[1] + dv[5]) / 2.0
dv3 = dv[2]
dv4 = dv[4]
pilist = [pi1, pi2, pi3, pi4]
dvlist = [dv1, dv2, dv3, dv4]
for i in range(4):
    plt.plot(pred, pilist[i]/dvlist[i], marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
xtick = np.arange(0.5, 1.01, 0.1)
ytick = np.arange(0.5, 1.41, 0.1)
plt.xticks(xtick, color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
plt.yticks(ytick)
plt.ylim(0.48, 1.42)
plt.xlim(0.48, 1.02)
plt.legend(loc='upper left')

plt.subplot(212)
for i in range(4):
    plt.plot(pred, dvlist[i], marker='o', ms=5, lw=0, alpha=0.75,
             label=labs[i], color=clist[i])
xtick = np.arange(0.5, 1.01, 0.1)
plt.xticks(xtick)
plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
plt.xlim(0.48, 1.02)
plt.ylabel(r'$D/\bar{D}$', labelpad=3)
ytick = np.arange(0.9, 1.31, 0.1)
plt.yticks(ytick)
plt.ylim(0.82, 1.38)
plt.legend(loc='upper left')
plt.savefig(f_save, dpi=512)
plt.close()

#%% NEUTRAL SITES PER BIN
f_save = root_dir + '/result/final_files/sfigs/sites_count_ratio.png'
plt.figure(figsize=(4,3))
plt.subplots_adjust(top=0.97, right=0.99, bottom=0.15, left=0.16)

plt.plot(snp[:,6] / np.sum(res[:,6:], axis=1))
plt.xlabel('B bin')
plt.ylabel('neutral polymorphism / neutral algined sites')
plt.savefig(f_save, dpi=512)
plt.close()
#%% CMMB AND CONSERVATION COMBINED
# # load results
# # f_sort = fdir + 'predsort_more.txt'
# f_save = root_dir + '/result/final_files/sfigs/{}_cmmb_and_cons.png'.format(fldr)
#
# # bbin, at_rate, a2t_rate, gc_rate, all_rate, hm_rate = np.loadtxt(f_bbin).T
# # div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
# plt.figure(figsize=(4, 7))
# plt.subplots_adjust(top=0.95, right=0.99, hspace=0.1, bottom=0.07, left=0.16)
# plt.subplot(211)
# plt.title('Confounder variables')
# plt.plot(pred, np.log10(cm * 1e6), marker='o', ms=5, lw=0, alpha=0.75,
#              color='deepskyblue')
# xtick = np.arange(0.5, 1.01, 0.1)
# plt.xticks(xtick, color='white')
# plt.ylabel('recombination rate (log10 cM/Mb)', labelpad=3)
# plt.ylim(-1.2, 1.2)
# plt.xlim(0.48, 1.02)
#
# plt.subplot(212)
# plt.plot(pred, cn / 1000.0, marker='o', ms=5, lw=0, alpha=0.75,
#          color='purple')
# xtick = np.arange(0.5, 1.01, 0.1)
# plt.xticks(xtick)
# plt.ylabel('conservation score', labelpad=3)
# plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
# plt.xlim(0.48, 1.02)
# plt.savefig(f_save, dpi=512)
# plt.close()

#%%