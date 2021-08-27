__author__ = 'davidmurphy'

import os
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir
from figures.sub_rate_figs import neutcons_uratio_single
from skmisc import loess
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from figures.common_functions import format_panels, cst_from_fldr


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'
figdir = final_dir + '/mainfigs'


# FUNCTIONS USED
def rsq_from_fldr(fldr):
    f_rsq = root_dir + '/result/final_files/{}/rsq.log'.format(fldr)
    r = np.loadtxt(f_rsq)
    return r


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def get_loess_line(fldr, span, return_points=False, load_con=False):
    # load results in 2000 bins for LOESS plots
    fdir = final_dir + '/{}/'.format(fldr)
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 2000)
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


#%% FIGURE 2B
def compare_rsq(flist, llist, clist, sname, sfldr, fletter='B'):
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(flist[0]))[:16,0])
    print w
    for fl in flist:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    plt.figure(figsize=(2.5,2.5))
    plt.subplots_adjust(right=1, top=1, left=0.25, bottom=0.17)
    for i in xrange(len(flist)):
        rs = rsq_list[i]
        lab = llist[i]
        col = clist[i]
        plt.plot(w, rs, label=lab, color=col, marker='o', lw=0, ms=6,
                 alpha=0.8)

    plt.xlabel('window size (log-scale)')
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4, 5, 6]])
    plt.ylim(0, 0.65)
    plt.ylabel('proportion of variance\n' + r'explained $(R^2)$')
    plt.legend(prop=dict(size=9), loc='upper left', handletextpad=0.1,
               borderaxespad=0.01)
    plt.text(0.01, 0.95, fletter, transform=plt.gcf().transFigure)

    f_save = final_dir + '/{}/{}.rsq.png'.format(sfldr, sname)
    plt.savefig(f_save, dpi=512)
    plt.close()

# Figure 2B
flist = ['ape_cons94_clean', 'ape_cons94_exonic']
clist = ['darkorange', 'fuchsia']
llist = ['all conserved', 'exonic conserved']
sname = 'fig_2B'
sfldr = 'mainfigs'
compare_rsq(flist, llist, clist, sname, sfldr, 'B')

#%% FIGURE 2A
def figure_2A(fldr):
    """Figure 2A plot"""
    # fixed variables:
    color = 'darkorange'
    label = 'predicted'

    # set chr1 file name
    fchr1 = 'chr1.1.00e+06win_0.5slide_pred_and_obs.txt'
    f_name = final_dir + '/{}/{}'.format(fldr, fchr1)

    # load and normalize data
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    # remove some outliers
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan

    # create new figure
    plt.figure(figsize=(6.4, 2))
    plt.subplots_adjust(left=0.07, bottom=0.2, right=1, top=0.9)

    # plot observed and predicted diversity
    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(xi, prd, color=color, lw=1.1, alpha=0.8, label=label)

    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    plt.ylim(0, 1.9)
    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    plt.xticks(xrange(25, 250, 50), y=0.05)
    plt.xlim(0, xi[-1])

    plt.legend(loc='lower center', ncol=2, frameon=1,
               framealpha=0.75, facecolor='white', handlelength=1.2)

    # figure letter
    plt.text(0.01, 0.93, 'A', transform=plt.gcf().transFigure)
    f_save = figdir + '/fig2A.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


# Figure 2A
fldr = 'cadd93'
figure_2A(fldr)


#%% FIGURE 2A-B
def figure_2AB(fldr, changebackground=False):
    """Figure 2A-B plot"""
    # fixed variables:
    color = 'darkorange'
    label = 'predicted'
    # fldr_1 = 'ape_cons94_clean'
    fldr_1 = fldr
    fldr_2 = 'fish_cons94_gmask_exonic'
    clist = ['darkorange', 'dodgerblue']
    llist = ['all', 'exonic']
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # set chr1 file name
    fchr1 = 'chr1.1.00e+06win_0.5slide_pred_and_obs.txt'
    f_name = final_dir + '/{}/{}'.format(fldr, fchr1)

    # load and normalize data
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    # remove some outliers
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan

    # load R^2 data for all conserved and exonic conserved (ape)
    rsq_list = []
    flist = [fldr_1, fldr_2]
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(fldr_1))[:16, 0])
    for fl in flist:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    # create new figure
    plt.figure(figsize=(6.4, 2))
    # plt.subplots_adjust(left=0.07, bottom=0.18, right=0.995, top=0.9, wspace=0.32)
    plt.subplots_adjust(left=0.07, bottom=0.17, right=0.998, top=0.99, wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1,2))
    format_panels(ax1)

    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(xi, prd, color=color, lw=1.1, alpha=0.8, label=label)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(0, 1.9)
    plt.ylim(0, 1.99)

    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x/25)%2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])

    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)
    format_panels(ax2)

    for i in xrange(len(flist)):
        rs = rsq_list[i]
        lab = llist[i]
        col = clist[i]
        plt.plot(w, rs, label=lab, color=col, marker='o', lw=0, ms=4, alpha=0.8)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $\mathrm{(R^2)}$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1,
               labelspacing=0.25)
    # figure letters
    # plt.text(0.36, 0.93, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.86, 0.93, 'B', transform=plt.gcf().transFigure)
    plt.text(0.075, 0.915, 'A', transform=plt.gcf().transFigure)
    plt.text(0.7475, 0.915, 'B', transform=plt.gcf().transFigure)

    if changebackground:
        f_save = figdir + '/fig2AB.{}.chr1.white.png'.format(fldr)
    else:
        f_save = figdir + '/fig2AB.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'fish_cons94_gmask_mnb_378'
figure_2AB(fldr, changebackground=True)


#%% FIGURE 2A-B FOR GUY
def figure_2AB_for_guy(fldr, changebackground=False):
    """Figure 2A-B plot"""
    # fixed variables:
    color = 'darkorange'
    label = 'predicted'
    # fldr_1 = 'ape_cons94_clean'
    fldr_1 = fldr
    fldr_2 = 'fish_cons94_gmask_exonic'
    clist = ['darkorange', 'dodgerblue']
    llist = ['all', 'exonic']
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # set chr1 file name
    fchr1 = 'chr1.1.00e+06win_0.5slide_pred_and_obs.txt'
    f_name = final_dir + '/{}/{}'.format(fldr, fchr1)

    # load and normalize data
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    # remove some outliers
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan

    # load R^2 data for all conserved and exonic conserved (ape)
    rsq_list = []
    flist = [fldr_1]
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(fldr_1))[:16, 0])
    for fl in flist:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    # create new figure
    plt.figure(figsize=(6.4, 2))
    # plt.subplots_adjust(left=0.07, bottom=0.18, right=0.995, top=0.9, wspace=0.32)
    plt.subplots_adjust(left=0.07, bottom=0.17, right=0.998, top=0.99, wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1,2))

    # version with white background
    if changebackground:
        ax1.set_facecolor('white')
        ax1.grid(**gridargs)
        ax1.axis('on')
        ax1.patch.set_edgecolor('black')
        ax1.patch.set_linewidth(axis_lw)

    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(xi, prd, color=color, lw=1.1, alpha=0.8, label=label)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(0, 1.9)
    plt.ylim(0, 1.99)

    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x/25)%2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])

    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)

    # version with white background
    if changebackground:
        ax2.set_facecolor('white')
        ax2.grid(**gridargs)
        ax2.axis('on')
        ax2.patch.set_edgecolor('black')
        ax2.patch.set_linewidth(axis_lw)

    for i in xrange(len(flist)):
        rs = rsq_list[i]
        lab = llist[i]
        col = clist[i]
        plt.plot(w, rs, label=lab, color=col, marker='o', lw=0, ms=4, alpha=0.8)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $\mathrm{(R^2)}$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    # plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
    #            borderaxespad=0.3, ncol=1, frameon=1,
    #            framealpha=0.75, facecolor='white', columnspacing=0.1,
    #            labelspacing=0.25)
    # figure letters
    # plt.text(0.36, 0.93, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.86, 0.93, 'B', transform=plt.gcf().transFigure)
    plt.text(0.075, 0.915, 'A', transform=plt.gcf().transFigure)
    plt.text(0.7475, 0.915, 'B', transform=plt.gcf().transFigure)

    if changebackground:
        f_save = figdir + '/fig2AB.{}.forGuy.chr1.white.png'.format(fldr)
    else:
        f_save = figdir + '/fig2AB.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_mnb_378'
figure_2AB_for_guy(fldr, changebackground=True)


#%% FIGURE 2A-B FINAL FOR PAPER
def figure_2AB_final_format(fldr):
    """Figure 2A-B plot"""
    # fixed variables:
    color = 'darkorange'
    label = 'predicted (CADD)'
    clist = ['darkorange', 'dodgerblue']

    # set chr1 file name
    fchr1 = 'chr1.1.00e+06win_0.5slide_pred_and_obs.txt'
    f_name = final_dir + '/{}/{}'.format(fldr, fchr1)

    # load and normalize data
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    # remove some outliers
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan

    # create new figure
    plt.figure(figsize=(6.4, 2))
    plt.subplots_adjust(left=0.07, bottom=0.17, right=0.998, top=0.99,
                        wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1,2))
    format_panels(ax1)

    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(xi, prd, color=color, lw=1.1, alpha=0.8, label=label)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(0, 1.9)
    plt.ylim(0, 1.99)

    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x/25)%2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])

    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)
    format_panels(ax2)

    # load R^2 data for all conserved and exonic conserved (ape)
    rsq_folders = [fldr, 'fish_cons94_new',
                   'cadd94_gmask_exonic','fish_cons94_new_exonic']
    rsq_shapes = ['o', '^', 'o', '^']
    rsq_labels = ['CADD', 'phastCons', r'$\mathrm{CADD_{e}}$',
                  r'$\mathrm{phastCons_{e}}$']
    rsq_colors = ['darkorange', 'deepskyblue', 'orangered', 'blue']
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(fldr))[:16, 0])
    for fl in rsq_folders:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    for i in xrange(len(rsq_folders)):
        rs = rsq_list[i]
        lab = rsq_labels[i]
        shape = rsq_shapes[i]
        col = rsq_colors[i]
        alpha = 0.6 if i in [1, 3] else 1.0
        plt.plot(w, rs, label=lab, color=col, marker=shape, lw=0, ms=5,
                 alpha=alpha)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1,
               labelspacing=0.25, prop=dict(size=8.5))

    plt.text(0.075, 0.915, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.7475, 0.915, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/updateMarch2021.fig2AB.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat'
figure_2AB_final_format(fldr)


#%% FIGURE 3
def figure_3(fldr, fcl, xlab, letter=None):
    # fixed variables:
    color = 'darkorange'
    label = 'predicted'

    # COLLATED PLOT
    f_col = final_dir + '/{}/{}.collate.5.00e-03width.npy'.format(fldr, fcl)
    bins, div, pi, pr, cnts = np.load(f_col).T
    obs = pi / cnts
    prd = pr / cnts
    y_max = obs.mean()
    obs /= y_max
    prd /= y_max
    bins += 2.5e-3

    plt.figure(figsize=(3.25, 1.95))
    plt.subplots_adjust(left=0.135, bottom=0.175, right=0.995, top=0.995)
    ax1 = plt.subplot(111)
    format_panels(ax1)

    plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(bins, prd, label=label, color=color, alpha=0.8, lw=1.1)
    plt.xlabel(xlab, labelpad=2)
    plt.xticks(y=0.05)
    # plt.xlim(-0.41, 0.41)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    if fcl == 'YRI_syn_s1':
        plt.ylim(0.81, 1.13)
    if fcl == 'hc_derived_cons':
        plt.yticks([0.7, 1.0, 1.3], x=0.03)
    else:
        plt.yticks(x=0.03)
    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)
    if letter is not None:
        plt.text(0.145, 0.915, letter, transform=plt.gcf().transFigure)

    f_save = figdir + '/fig3.{}.{}.collated.png'.format(fldr, fcl)
    plt.savefig(f_save, dpi=512)
    plt.close()


def figure_3_hist(fldr, fcl, xlab, changebackground):
    # fixed variables:
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # COLLATED PLOT
    f_col = final_dir + '/{}/{}.collate.5.00e-03width.npy'.format(fldr, fcl)
    bins, div, pi, pr, cnts = np.load(f_col).T
    bins += 2.5e-3

    plt.figure(figsize=(5, 3))
    plt.subplots_adjust(left=0.12, bottom=0.15, right=1, top=1)
    ax1 = plt.subplot(111)

    # version with white background
    if changebackground:
        ax1.set_facecolor('white')
        ax1.grid(**gridargs)
        ax1.axis('on')
        ax1.patch.set_edgecolor('black')
        ax1.patch.set_linewidth(axis_lw)

    plt.plot(bins, np.log10(cnts), color='k', lw=0.8)
    plt.xlabel(xlab, labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'neutral site counts (log10 scale)', labelpad=2)
    plt.yticks(x=0.02)
    # plt.legend(loc='lower right', ncol=2, frameon=1,borderaxespad=0.3,
    #            borderpad=0.2, framealpha=0.75, facecolor='white',
    #            handlelength=1.2)

    if changebackground:
        f_save = figdir + '/fig3.{}.{}.counts.white.png'.format(fldr, fcl)
    else:
        f_save = figdir + '/fig3.{}.{}.counts.png'.format(fldr, fcl)
    plt.savefig(f_save, dpi=512)
    plt.close()


# focals = ['YRI_syn_s1', 'hc_derived_cons', 'exon', 'cadd94_gmask_exonic']
# xlabels = ['synonymous\nsubstitution', 'conserved\nsubstitution', 'exon',
#            'conserved\nexonic segment']
# letters = 'A B C D'.split()
# idx = 0
# indices = range(4)
# for idx in indices:
#     xlab = 'distance to nearest {} (cM)'.format(xlabels[idx])
#     fcl = focals[idx]
#     figure_3(fldr, fcl, xlab, letters[idx])

fldr = 'cadd94_gmask_v1.6_without_bstat'
xlab = 'distance to nearest NS substitution (cM)'
figure_3(fldr, 'nonsyn', xlab)



#%% FIGURE 4 DATA
def figure_4_data():
    pct = [95, 94, 93]
    s_ape = [1 - neutcons_uratio_single('ape', pc=p)[0].mean() for p in pct]
    s_cad = [1 - neutcons_uratio_single('cadd', pc=p)[0].mean() for p in pct]
    s_fis = [1 - neutcons_uratio_single('fish', pc=p)[0].mean() for p in pct]

    # get inferred udel for ape and CADD
    u_fis = []
    u_cad = []
    u_ape = []
    for p in pct:
        if p == 94:
            cadfldr = 'cadd94_gmask_v1.6_without_bstat'
            fisfldr = 'fish_cons94_new'
            apefldr = 'ape_cons94_gmask_mnb_378'

        else:
            cadfldr = 'cadd{}_gmask_v1.6_without_bstat'.format(p)
            # cadfldr = 'cadd{}_gmask'.format(p)
            fisfldr = 'fish_cons{}_new'.format(p)
            apefldr = 'ape_cons{}_gmask'.format(p)

        u_cad.append(cst_from_fldr(cadfldr).stat.utot[0] * 1e8)
        u_fis.append(cst_from_fldr(fisfldr).stat.utot[0] * 1e8)
        u_ape.append(cst_from_fldr(apefldr).stat.utot[0] * 1e8)

    # return s_ape, s_cad, u_ape, u_cad
    return s_ape, s_cad, s_fis, u_ape, u_cad, u_fis


s_ape, s_cad, s_fis, u_ape, u_cad, u_fis = figure_4_data()


#%% FIGURE 4
def figure_4(s_fis, s_cad, u_fis, u_cad, changebackground):
    #  upper and lower bounds of u total
    ulo = 1.29
    uhi = 1.51
    pct = [95, 94, 93]

    # axis format variables for white background
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # plt.figure(figsize=(5, 3))
    fig = plt.figure(figsize=(3.25, 1.95))
    plt.subplots_adjust(wspace=0.05, left=0.2, right=1, top=0.995, bottom=0.21)
    # panel A: ape rates
    ax1 = plt.subplot(121)
    # version with white background
    if changebackground:
        ax1.set_facecolor('white')
        ax1.grid(**gridargs)
        ax1.axis('on')
        ax1.patch.set_edgecolor('black')
        ax1.patch.set_linewidth(axis_lw)

    x1 = np.arange(len(pct))
    utop = np.array([u / ulo for u in u_fis])
    ubot = np.array([u / uhi for u in u_fis])
    uheight = utop-ubot
    # print utop
    # print ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_fis, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('phastCons', labelpad=2)
    plt.ylabel('proportion of\ndeleterious mutations', labelpad=2)
    plt.yticks(x=0.06)
    plt.ylim(0, 1.15)

    # panel B: CADD rates
    ax2 = plt.subplot(122)
    # version with white background
    if changebackground:
        ax2.set_facecolor('white')
        ax2.grid(**gridargs)
        ax2.axis('on')
        ax2.patch.set_edgecolor('black')
        ax2.patch.set_linewidth(axis_lw)

    utop = np.array([u / ulo for u in u_cad])
    ubot = np.array([u / uhi for u in u_cad])
    uheight = utop-ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_cad, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('CADD', labelpad=2)
    plt.yticks(color='none', x=0.05)
    plt.ylim(0, 1.15)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=6.5))
    # leg._legend_box.align = 'left'
    # plt.setp(leg.get_title(), fontsize='xx-small')
    # customize legend
    handles, labels = ax1.get_legend_handles_labels()
    leg = fig.legend(reversed(handles), reversed(labels), loc='lower center',
                     title='estimates from:', frameon=1, framealpha=1,
                     facecolor='white', prop=dict(size=8.5), ncol=2,
                     borderpad=0.2, labelspacing=0.2, handlelength=1.1,
                     handletextpad=0.4, bbox_to_anchor=(0.6, 0.195),
                     columnspacing=0.8)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=9),
    #                  ncol=2)
    leg._legend_box.align = 'left'
    leg.set_title('estimates from:', prop={'size':8.5})
    plt.text(0.21, 0.92, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.62, 0.92, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/updateMarch2021.fig4.urates.png'

    plt.savefig(f_save, dpi=512)
    plt.close()


figure_4(s_fis, s_cad, u_fis, u_cad, changebackground=True)


#%% FIGURE 4 FOR GUY
def figure_4_for_guy(s_fis, s_cad, u_fis, u_cad, changebackground):
    #  upper and lower bounds of u total
    ulo = 1.29
    uhi = 1.51
    pct = [95, 94, 93]
    x1 = np.arange(len(pct))

    # axis format variables for white background
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    fig = plt.figure(figsize=(2.5, 2.5))
    plt.subplots_adjust(wspace=0.05, left=0.27, right=1, top=0.995, bottom=0.4)

    # panel B: CADD rates
    ax2 = plt.subplot(111)
    # version with white background
    if changebackground:
        ax2.set_facecolor('white')
        ax2.grid(**gridargs)
        ax2.axis('on')
        ax2.patch.set_edgecolor('black')
        ax2.patch.set_linewidth(axis_lw)

    utop = np.array([u / ulo for u in u_cad])
    ubot = np.array([u / uhi for u in u_cad])
    uheight = utop-ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_cad, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('CADD', labelpad=2)
    plt.ylabel('proportion of\ndeleterious mutations', labelpad=2)
    plt.yticks(x=0.04)
    plt.ylim(0, 1.15)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=9))
    handles, labels = ax2.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc='lower center',
                     title='estimates from:', frameon=1, framealpha=1,
                     facecolor='white', prop=dict(size=10), ncol=1,
                     borderpad=0.2, labelspacing=0.2, handlelength=1.1,
                     handletextpad=0.4, bbox_to_anchor=(0.6, -0.01),
                     columnspacing=0.8)

    if changebackground:
        f_save = figdir + '/fig4.urates.white.forGuy.png'
    else:
        f_save = figdir + '/fig4.urates.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


figure_4_for_guy(s_fis, s_cad, u_fis, u_cad, changebackground=True)


#%% FIGURE 5
def figure_5(fldr, changebackground):
    # fixed variables
    num = 100
    axmin, axmax = 0.55, 1.12
    color = 'darkorange'
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # load data
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, num)
    div, pi, pred = np.loadtxt(sort_file).T
    rst = cst_from_fldr(fldr)
    # normalize by pi0
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)

    # version with white background
    if changebackground:
        ax1.set_facecolor('white')
        ax1.grid(**gridargs)
        ax1.axis('on')
        ax1.patch.set_edgecolor('black')
        ax1.patch.set_linewidth(axis_lw)

    # plot predicted vs. observed
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor=color, markeredgewidth=0.9, lw=0,
             alpha=0.75)
    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='darkslategray', ls='--', alpha=0.65)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(axmin+0.02, 1.02, 'without linked selection', ha='left',
             va='center')

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    xtick = np.arange(0.5, 1.2, 0.1)
    plt.xticks(xtick, y=0.02)
    plt.ylim(axmin, axmax)
    plt.xlim(axmin, axmax)
    # plt.legend(loc='lower right', ncol=1, frameon=1,
    #            framealpha=0.75, facecolor='white')
    plt.legend(loc='lower right')

    if changebackground:
        f_save = figdir + '/fig5.{}.sort.white.png'.format(fldr)
    else:
        f_save = figdir + '/fig5.{}.sort.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


figure_5('cadd94_gmask_mnb_378', True)


#%% FIGURE 5 DETAIL
def figure_5_detail(fldr, pct, span, changebackground):
    # fixed variables
    num = 2000
    # axmin, axmax = 0.55, 1.12
    color = 'darkorange'
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # load data
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, num)
    div, pi, pred = np.loadtxt(sort_file).T
    rst = cst_from_fldr(fldr)
    # normalize by pi0
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # create LOESS line to put alongside detail region
    weights = np.ones(shape=len(pi))
    # span = 0.02
    pr_loess = predict_loess(pred, pi, weights, span, pred)

    # get the index corresponding to the percentile of sites to take
    # idx = int(pct * num)
    idx = np.where(pred>=pct)[0][0]
    print 1.0 - idx / 2000.0
    # take only the top percentile of sites
    pi = pi[idx:]
    pred = pred[idx:]
    pr_loess = pr_loess[idx:]

    pmin = pred.min()
    pmax = pred.max()
    axspan = pmax - pmin
    axmin = pmin-0.05*axspan
    axmax = pmax+0.05*axspan

    # create new plot
    plt.figure(figsize=(2.5, 2.5))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.18, bottom=0.16)
    ax1 = plt.subplot(111)

    # version with white background
    if changebackground:
        ax1.set_facecolor('white')
        ax1.grid(**gridargs)
        ax1.axis('on')
        ax1.patch.set_edgecolor('black')
        ax1.patch.set_linewidth(axis_lw)

    # plot predicted vs. observed
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor=color, markeredgewidth=0.9, lw=0,
             alpha=0.75)
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='darkslategray', ls='--', alpha=0.65)
    loess_lab = 'LOESS fit'
    # loess_lab = 'LOESS fit\n(span={})'.format(str(span)
    plt.plot(pred, pr_loess, label=loess_lab, color='mediumvioletred',
             alpha=0.65)
    # plt.text(0.52, 1.17, r'$n_{bins}=$' + str(num))
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.yticks(x=0.02)
    plt.ylim(0.79, 1.29)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    plt.xlim(axmin, axmax)
    plt.xticks(y=0.02)
    plt.legend(loc='upper left')

    if changebackground:
        f_save = figdir + '/fig5_detail.{}.sort.white.png'.format(fldr)
    else:
        f_save = figdir + '/fig5_detail.{}.sort.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


figure_5_detail('cadd93_bth600', 0.9, 0.01, True)


#%% FIGURE 5 WITH INSET
def figure_5_inset(fldr, pct, span):
    # fixed variables
    num = 100
    axmin, axmax = 0.55, 1.12
    color = 'darkorange'

    # load data
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, num)
    div, pi, pred = np.loadtxt(sort_file).T
    rst = cst_from_fldr(fldr)
    # normalize by pi0
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # INSET DATA
    inum = 2000
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, inum)
    _, ipi, ipred = np.loadtxt(sort_file).T
    # normalize by pi0
    ipi /= pi0
    ipred /= pi0
    # get the index corresponding to the percentile of sites to take
    # idx = int(pct * num)
    idx = np.where(ipred>=pct)[0][0]
    print 1.0 - idx / 2000.0
    # take only the top percentile of sites
    ipi = ipi[idx:]
    ipred = ipred[idx:]
    # pr_loess = pr_loess[idx:]

    # create LOESS line to put alongside detail region
    weights = np.ones(shape=len(ipi))
    # span = 0.02
    pr_loess = predict_loess(ipred, ipi, weights, span, ipred)

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)

    # plot predicted vs. observed
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor=color, markeredgewidth=0.9, lw=0,
             alpha=0.75)
    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='darkslategray', ls='--', alpha=0.65)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(axmin+0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    # BOX ON INSET DATA
    inset_lw = 0.5
    rect = Rectangle(xy=(pct, pct), width=0.1, height=0.215,
                     linewidth=inset_lw, edgecolor='red',
                     facecolor='none')
    ax1.add_patch(rect)

    # # INSET CONNECTOR LINES
    # # xst, xen = 0.8735, 1.1132
    # xst, xen = 0.88, 1.1
    # xpts = [xst, 0.9]
    # ypts = [0.796, 0.9]
    # plt.plot(xpts, ypts, color='k', ls='-', lw=0.25)
    # xpts = [1, xen]
    # ypts = [0.9, 0.796]
    # plt.plot(xpts, ypts, color='k', ls='-', lw=0.25)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    xtick = np.arange(0.5, 1.2, 0.1)
    plt.xticks(xtick, y=0.02)
    plt.ylim(axmin, axmax)
    plt.xlim(axmin, axmax)
    # plt.legend(loc='lower center', bbox_to_anchor=(0.3, 0.05))
    plt.text(0.75, 0.81, r'$y=x$', rotation=45, ha='center', va='center',
             color='darkslategray', alpha=0.65)

    # INSET FIGURE
    loess_lab = 'LOESS fit'
    # axins = inset_axes(ax1, width='42%', height='42%', loc=4,
    #                    borderpad=0.25)
    axins = inset_axes(ax1, width="100%", height="100%", loc=4,
                   bbox_to_anchor=(0.58, 0.03, 0.4, 0.4),
                       bbox_transform=ax1.transAxes)
    axins.plot(ipred, ipi, marker='o', ms=3.5, markerfacecolor='None',
             markeredgecolor=color, markeredgewidth=0.5, lw=0,
             alpha=0.5)
    axins.plot(ipred, ipred, color='darkslategray', ls='--', alpha=0.65)
    axins.axhline(y=1, color='k', alpha=0.8, ls='-')
    axins.plot(ipred, pr_loess, label=loess_lab, color='mediumvioletred',
             alpha=0.65)
    axins.legend(prop=dict(size=8))
    # axins.set_title('1Mb detail')
    # axins.set_xticks([0.9, 1.0], color='red', size=8)
    plt.xticks([0.9, 1.0], color='r', size=8, y=0.07)
    plt.xlim(0.9, 1.0)
    axins.set_yticks([])
    axins.set_facecolor('white')
    axins.grid(lw=0)
    axins.axis('on')
    axins.patch.set_edgecolor('red')
    axins.patch.set_linewidth(inset_lw)

    f_save = figdir + '/fig5.{}.sort.inset.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


figure_5_inset('cadd94_gmask_mnb_378', 0.9, 0.1)


#%% OBSERVED VS. PREDICTED FINAL FORMATING
def figure_5_final_format(fldr, span=0.1):
    """create loess smoothed plots for conservation and recombination rates"""
    # get 100 points data
    fdir = final_dir + '/{}/'.format(fldr)
    # f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
    # div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 100)
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
# fldr = 'cadd94_gmask_v1.6'

figure_5_final_format(fldr)


#%% standard figure set
fldr = 'fish_cons94_gmask_YRI_LD'
figure_2AB_for_guy(fldr)
figure_3(fldr, 'nonsyn', 'cM to nearest AA substitution')
figure_5_inset(fldr, 0.9, 0.1)


#%%
for gm in ['YRI_LD', 'deCODE_2019']:
    for an in ['fish_cons94', 'cadd94']:
        fldr = '{}_gmask_{}'.format(an, gm)
        figure_3(fldr, 'nonsyn', 'cM to nearest AA substitution')


#%%