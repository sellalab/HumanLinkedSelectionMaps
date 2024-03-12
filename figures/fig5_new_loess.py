__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir
from loess import loess_1d
from figures.common_functions import format_panels, cst_from_fldr


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = r'/Users/MURPHYD/Dropbox (OMRF)/final_files'
# final_dir = root_dir + '/result/final_files'
figdir = final_dir + '/mainfigs'


# FUNCTIONS USED
def rsq_from_fldr(fldr):
    f_rsq = root_dir + '/result/final_files/{}/rsq.log'.format(fldr)
    r = np.loadtxt(f_rsq)
    return r


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    # lo = loess.loess(xi, yi, weights=wts, span=span)
    lo = loess_1d.loess_1d(xi, yi, xnew=xtest, frac=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def get_loess_line(fldr, span, return_points=False, load_con=False):
    # load results in 2000 bins for LOESS plots
    fdir = final_dir + '/{}/'.format(fldr)
    # sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 2000)
    sort_file = final_dir + '/{}/basic_sort_n{}_leaveOneOut.txt'.format(fldr, 2000)
    div, pi, pred = np.loadtxt(sort_file).T
    # f_sort = fdir + 'sort_gc_cm_cn_il_n2000.txt'
    # div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    # load sorted conserved separately (done using 0.05cM radius of sites)
    if load_con:
        fcon = final_dir + '/{}/sorted.cons.n2000.txt'.format(fldr)
        cn = np.loadtxt(fcon)[:2000,2]

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # points = [pi, gc, cn, cm, il]
    points = [pi]
    wts = np.ones(shape=len(pred))
    loess_lines = []
    for a in points:
        lo_line = predict_loess(pred, a, wts, span, pred)
        loess_lines.append(lo_line)

    if return_points:
        return pred, loess_lines, points
    else:
        return pred, loess_lines


#%% OBSERVED VS. PREDICTED FINAL FORMATING
def figure_5_final_format(fldr, span=0.1):
    """create loess smoothed plots for conservation and recombination rates"""
    # get 100 points data
    fdir = final_dir + '/{}/'.format(fldr)
    # f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
    # div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T
    # sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 100)
    sort_file = final_dir + '/{}/basic_sort_n{}_leaveOneOut.txt'.format(fldr, 2000)
    div, pi, pred = np.loadtxt(sort_file).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # get loess line and original points
    prlo, lolist = get_loess_line(fldr, span)
    # pilo, gclo, cnlo, cmlo, illo = lolist
    pilo = lolist[0]

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)

    axmin, axmax = 0.55, 1.2
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.94, 0.48
    xtick = np.arange(0.5, 1.2, 0.1)

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='k', ls='--')
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    # plot predicted vs. observed
    plt.plot(pred, pi, marker='o', ms=5, color='darkorange', lw=0,
             alpha=0.5)
    # plot LOESS line
    plt.plot(prlo, pilo, lw=2, color='orangered')
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xticks(xtick, y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.52, axmax)
    plt.xlim(axmin, 1.02)
    # solve y=x rotation
    adlen = axmax - 0.52
    oplen = adlen * (1.02 - axmin) / (axmax - 0.52)
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')
    # plt.text(0.75, 0.81, r'$y=x$', rotation=38, ha='center', va='center',
    #          color='darkslategray', alpha=0.65)
    # plt.text(tx1, ty1, 'A', transform=plt.gcf().transFigure)

    f_save = figdir + '/updateMarch2021.fig5.{}.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


# fldr = 'cadd94_gmask_mnb_378'
fldr = 'cadd94_gmask_v1.6_without_bstat'
# fldr = 'cadd94_gmask_v1.6_without_bstat_jackknife_results'
# fldr = 'cadd94_gmask_v1.6'

figure_5_final_format(fldr)

#%%
jkconstpred = np.load(final_dir+'/cadd94_gmask_v1.6_without_bstat_jackknife_results/cadd94_gmask_v1.6_without_bstat.jkpred_const.npy')
jkpred = np.load(final_dir+'/cadd94_gmask_v1.6_without_bstat_jackknife_results/cadd94_gmask_v1.6_without_bstat.jkpred.npy')

#%%