__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir
from figures.sub_rate_figs import neutcons_uratio_single
from skmisc import loess
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from figures.common_functions import format_panels, cst_from_fldr

# set the seaborn pallete
# seaborn.set_style()
# seaborn.set_theme({'axes.facecolor': 'white', 'grid.color': 'k', 'grid.linestyle': '--', 'grid.linewidth': 0.2,
#                    'axes.edgecolor': 'k'})
# seaborn.set()
style_dict = {'axes.facecolor': 'white', 'grid.color': 'k', 'grid.linestyle': '--', 'grid.linewidth': 0.2,
              'axes.edgecolor': 'k', 'axes.linewidth': 0.25, 'legend.fancybox': True, 'grid.alpha': 0.5,
              'lines.dashed_pattern': [10, 10], 'lines.scale_dashes': 1}
seaborn.set_theme(style_dict)
for (k, v) in style_dict.items():
    plt.rcParams[k] = v

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = r'/Users/MURPHYD/Dropbox (OMRF)/linked_selection/result/final_files'
# final_dir = root_dir + '/result/final_files'
figdir = final_dir + '/mainfigs'


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def get_loess_line(fldr, span, return_points=False, load_con=False, loo=False):
    # load results in 2000 bins for LOESS plots
    fdir = final_dir + '/{}/'.format(fldr)
    if loo:
        sort_file = final_dir + '/{}/basic_sort_n{}_leaveOneOut.txt'.format(fldr, 2000)
    else:
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 2000)
    div, pi, pred = np.loadtxt(sort_file).T
    # f_sort = fdir + 'sort_gc_cm_cn_il_n2000.txt'
    # div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    # load sorted conserved separately (done using 0.05cM radius of sites)
    if load_con:
        fcon = final_dir + '/{}/sorted.cons.n2000.txt'.format(fldr)
        cn = np.loadtxt(fcon)[:2000, 2]

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


#%% FIGURE 2A-B FINAL FOR PAPER
def figure_2AB_final_format(fldr):
    """Figure 2A-B plot"""
    # style_dict = {'axes.facecolor': 'white', 'grid.color': 'k', 'grid.linestyle': '--', 'grid.linewidth': 0.2,
    #               'axes.edgecolor': 'k', 'axes.linewidth': 0.25, 'legend.fancybox': True, 'grid.alpha': 0.5,
    #               'lines.dashed_pattern': [10, 10], 'lines.scale_dashes': 1}
    # seaborn.set_theme(style_dict)
    # for (k, v) in style_dict.items():
    #     plt.rcParams[k] = v

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
    plt.subplots_adjust(left=0.075, bottom=0.185, right=0.998, top=0.99,
                        wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1, 2))

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
    # plt.tick_params(axis=u'both', which=u'both', length=0)

    legend = plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.2,
                        borderpad=0.2, framealpha=0.8, facecolor='white',
                        handlelength=1.2, labelspacing=0.25)
    legend.get_frame().set_linewidth(0.5)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)

    # load R^2 data for all conserved and exonic conserved (ape)
    rsq_folders = [fldr, 'fish_cons94_new_jackknife_results',
                   'cadd94_gmask_exonic', 'fish_cons94_new_exonic']
    rsq_shapes = ['o', '^', 'o', '^']
    rsq_labels = ['CADD', 'phastCons', r'$\mathrm{CADD_{e}}$',
                  r'$\mathrm{phastCons_{e}}$']
    rsq_colors = ['darkorange', 'deepskyblue', 'orangered', 'blue']
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(rsq_folders[-1]))[:16, 0])
    for fl in rsq_folders:
        if 'jackknife' in fl:
            fsq_file = final_dir + '/{}/{}.jkrsq.log'.format(fl, fl.replace('_jackknife_results', ''))
        else:
            fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    for i in range(len(rsq_folders)):
        rs = rsq_list[i]
        lab = rsq_labels[i]
        shape = rsq_shapes[i]
        col = rsq_colors[i]
        alpha = 0.6 if i in [1, 3] else 1.0
        plt.plot(w, rs, label=lab, color=col, marker=shape, lw=0, ms=5,
                 alpha=alpha)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.0f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    # plt.tick_params(axis=u'both', which=u'both', length=0)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    legend = plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
                        borderaxespad=0.2, ncol=1, frameon=1,
                        framealpha=0.8, facecolor='white', columnspacing=0.1,
                        labelspacing=0.25, prop=dict(size=8.5))
    legend.get_frame().set_linewidth(0.5)
    fig_letters = ['A', 'B']
    # fig_letters = ['a', 'b']

    plt.text(0.08, 0.915, fig_letters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.7475, 0.915, fig_letters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/NatGenet.fig2AB.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


# fldr = 'cadd94_gmask_v1.6_without_bstat_jackknife_results'
fldr = 'mockYRI'
# fldr = 'dropped_chrom_results'
figure_2AB_final_format(fldr)


#%% FIGURE 2AB SIMULATED
def figure_2AB_simulated_data(fldr, label, label2, bxpad):
    """Figure 2A-B plot"""
    # style_dict = {'axes.facecolor': 'white', 'grid.color': 'k', 'grid.linestyle': '--', 'grid.linewidth': 0.2,
    #               'axes.edgecolor': 'k', 'axes.linewidth': 0.25, 'legend.fancybox': True, 'grid.alpha': 0.5,
    #               'lines.dashed_pattern': [10, 10], 'lines.scale_dashes': 1}
    # seaborn.set_theme(style_dict)
    # for (k, v) in style_dict.items():
    #     plt.rcParams[k] = v

    # fixed variables:
    color = 'darkorange'

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
    plt.subplots_adjust(left=0.075, bottom=0.185, right=0.998, top=0.99,
                        wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1, 2))

    plt.plot(xi, obs, label='observed (simulated data)', color='darkslategray', lw=1.1)
    plt.plot(xi, prd, color=color, lw=1.1, alpha=0.8, label='predicted (' + label + ')')
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(0, 1.9)
    plt.ylim(0, 1.99)

    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x/25)%2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])
    # plt.tick_params(axis=u'both', which=u'both', length=0)

    legend = plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.2,
                        borderpad=0.2, framealpha=0.8, facecolor='white',
                        handlelength=1.2, labelspacing=0.25)
    legend.get_frame().set_linewidth(0.5)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)

    # load R^2 data for all conserved and exonic conserved (ape)
    rsq_folders = [fldr, 'cadd94_gmask_v1.6_without_bstat']
    rsq_shapes = ['o', '^', 'o', '^']
    rsq_labels = [label, label2]
    rsq_colors = ['darkorange', 'deepskyblue', 'orangered', 'blue']
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(rsq_folders[-1]))[:16, 0])
    for fl in rsq_folders:
        if 'jackknife' in fl:
            fsq_file = final_dir + '/{}/{}.jkrsq.log'.format(fl, fl.replace('_jackknife_results', ''))
        else:
            fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    for i in range(len(rsq_folders)):
        rs = rsq_list[i]
        lab = rsq_labels[i]
        shape = rsq_shapes[i]
        col = rsq_colors[i]
        alpha = 0.6 if i in [1, 3] else 1.0
        plt.plot(w, rs, label=lab, color=col, marker=shape, lw=0, ms=5,
                 alpha=alpha)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.0f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    # plt.tick_params(axis=u'both', which=u'both', length=0)
    plt.ylim(-0.01, 0.65)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    legend = plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
                        borderaxespad=bxpad, ncol=1, frameon=1,
                        framealpha=0.8, facecolor='white', columnspacing=0.1,
                        labelspacing=0.25, prop=dict(size=8.5))
    legend.get_frame().set_linewidth(0.5)
    fig_letters = ['A', 'B']
    # fig_letters = ['a', 'b']

    plt.text(0.08, 0.915, fig_letters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.7475, 0.915, fig_letters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/reply2rev.fig2AB.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


# fldr = 'cadd94_gmask_v1.6_without_bstat_jackknife_results'
fldr = 'mockYRI'
lbl = 'CADD: simulated data'
lbl2 = 'CADD: real data'
bxpad = 0.66
# fldr = 'dropped_chrom_results'
# lbl = 'CADD: leave chrom out'
# lbl2 = 'CADD: all data'
# bxpad = 0.3
figure_2AB_simulated_data(fldr, lbl, lbl2, bxpad)


#%% FIGURE 2A-B REPLY TO REVIEWERS
def figure_2AB_reply2reviews(flist, clist, llist, lwdths, rsq_folders, rsq_labels, sname):
    """Figure 2A-B plot"""
    # fixed variables:
    color = 'darkorange'

    # get data for all the folders
    xilist, prdlist, obslist = [], [], []
    for fldr in flist:
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
        xilist.append(xi)
        prdlist.append(prd)
        obslist.append(obs)

    # create new figure
    plt.figure(figsize=(6.4, 2))
    plt.subplots_adjust(left=0.075, bottom=0.185, right=0.998, top=0.99,
                        wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1, 2))
    plt.plot(xilist[0], obslist[0], label='observed', color='darkslategray', lw=1.1)
    for i in range(len(prdlist)):
        plt.plot(xilist[i], prdlist[i], color=clist[i], lw=lwdths[i], alpha=0.6, label=llist[i])
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    plt.ylim(0, 1.99)

    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x/25)%2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xilist[0][-1])
    legend = plt.legend(loc='lower right', ncol=2, frameon=1, borderaxespad=0.2,
                        borderpad=0.2, framealpha=0.8, facecolor='white',
                        handlelength=1.2, labelspacing=0.25)
    legend.get_frame().set_linewidth(0.5)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)

    # load R^2 data for all conserved and exonic conserved (ape)
    # rsq_folders = [fldr, 'fish_cons94_new_jackknife_results',
    #                'cadd94_gmask_exonic', 'fish_cons94_new_exonic']
    rsq_shapes = ['D', '^', 'o', '^']
    # rsq_labels = ['CADD', 'phastCons', r'$\mathrm{CADD_{e}}$',
    #               r'$\mathrm{phastCons_{e}}$']
    rsq_colors = clist
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(rsq_folders[-1]))[:16, 0])
    for fl in rsq_folders:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = abs(np.loadtxt(fsq_file)[:16, 1])
        rsq_list.append(rs)

    for i in range(len(rsq_folders)):
        rs = rsq_list[i]
        lab = rsq_labels[i]
        shape = rsq_shapes[i]
        col = rsq_colors[i]
        alpha = 0.6
        plt.plot(w, rs, label=lab, color=col, marker=shape, lw=0, ms=5,
                 alpha=alpha)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.0f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.ylim(-0.01, 0.65)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    legend = plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
                        borderaxespad=0.2, ncol=1, frameon=1,
                        framealpha=0.8, facecolor='white', columnspacing=0.1,
                        labelspacing=0.25, prop=dict(size=8.5))
    legend.get_frame().set_linewidth(0.5)
    # fig_letters = ['A', 'B']
    fig_letters = ['a', 'b']

    plt.text(0.08, 0.915, fig_letters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.7475, 0.915, fig_letters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/reply2review.fig2AB.{}.chr1.png'.format(sname)
    plt.savefig(f_save, dpi=512)
    plt.close()


# chr1_folder = 'dropped_chrom_results'
# chr1_lb = 'predicted (CADD: leave chrom out)'
# rsq_fl = ['cadd94_gmask_v1.6_without_bstat', 'dropped_chrom_results']
# rsq_lb = ['CADD: all data', 'CADD: leave chrom out']
# figure_2AB_reply2reviews(chr1_folder, chr1_lb, rsq_fl, rsq_lb)
flist = ['cadd94_gmask_v1.6_without_bstat', 'cadd94_gmask_v1.6_without_bstat_jackknife_results',
         'dropped_chrom_results']
llist = ['all data', 'leave 2 Mb out', 'leave 1 autosome out']
clist = ['darkorange', 'deepskyblue', 'fuchsia']
lwidths = [2.66, 1.33, 0.66]
sname = 'compare_LOO'
flist = ['cadd94_gmask_v1.6_without_bstat', 'mockYRI']
figure_2AB_reply2reviews(flist, clist, llist, lwidths, flist, llist, sname)


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
    plt.subplots_adjust(left=0.165, bottom=0.175, right=0.995, top=0.995)
    ax1 = plt.subplot(111)
    # format_panels(ax1)

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
    legend = plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
                        borderpad=0.2, framealpha=0.75, facecolor='white',
                        handlelength=1.2, labelspacing=0.25)
    legend.get_frame().set_linewidth(0.5)
    if letter is not None:
        plt.text(0.145, 0.915, letter, transform=plt.gcf().transFigure)

    f_save = figdir + '/fig3.{}.{}.collated.png'.format(fldr, fcl)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat_jackknife_results'
# fldr = 'mockYRI'
# fldr = 'dropped_chrom_results'
xlab = 'distance to nearest NS substitution (cM)'
figure_3(fldr, 'nonsyn', xlab)


#%% figure 5 special format for neutral data
def figure_5_neutral_coal_data(fldr, span=0.1):
    """create loess smoothed plots for conservation and recombination rates"""
    sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 100)
    div, pi, pred = np.loadtxt(sort_file).T
    # normalize by pi0
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0
    # get loess line and original points
    prlo, lolist = get_loess_line(fldr, span)
    pilo = lolist[0]
    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.19, bottom=0.13)
    ax1 = plt.subplot(111)
    axmin = pred.min()*0.99
    axmax = pred.max()*1.01
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
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(pi.min()*0.99, pi.max()*1.01)
    plt.xlim(pred.min()*0.999, 1.001)
    # solve y=x rotation
    adlen = axmax - 0.52
    oplen = adlen * (1.02 - axmin) / (axmax - 0.52)
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')
    f_save = figdir + '/fig5-neutral-coal-data.{}.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'mockYRI'
figure_5_neutral_coal_data(fldr)


#%% OBSERVED VS. PREDICTED FINAL FORMATING
def figure_5_final_format(fldr, span=0.1, loo=False):
    """create loess smoothed plots for conservation and recombination rates"""
    # get 100 points data
    plt.rcParams['hatch.linewidth'] = 0.75
    if loo:
        sort_file = final_dir + '/{}/basic_sort_n{}_leaveOneOut.txt'.format(fldr, 100)
    else:
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
    # format_panels(ax1)

    axmin, axmax = 0.55, 1.19
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.94, 0.48
    xtick = np.arange(0.5, 1.2, 0.1)

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='k', ls='--', dashes=(5, 3))
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    # plot 5 bins affected by thresholding as hollow
    plt.scatter(pred[:5], pi[:5], marker='o', s=25, edgecolors='darkorange', linestyle='-',
                hatch='//'*6, alpha=0.4, color='none', linewidths=0.75)
    # plt.plot(pred[:5], pi[:5], marker='o', ms=5, color='gray', lw=0,
    #          alpha=0.5, markerfacecolor='None')
    # plot remaining 95 bins as solid
    plt.plot(pred[5:], pi[5:], marker='o', ms=5, color='darkorange', lw=0,
             alpha=0.5)
    # make a big dot where 5 bins are averaged
    plt.scatter(np.mean(pred[:5]), np.mean(pi[:5]), s=30, alpha=0.9, marker='o', color='darkred', zorder=3)
    # print(pi[0:10])
    print(pred[:5])
    print(pred[-5:])
    # plot LOESS line
    plt.plot(prlo, pilo, lw=2, color='orangered')

    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xticks(xtick, y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.45, axmax)
    plt.xlim(axmin, 1.02)
    # plt.xlim(0.45, 1.45)
    # solve y=x rotation
    adlen = axmax - 0.52
    oplen = adlen * (1.02 - axmin) / (axmax - 0.52)
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')
    # plt.text(0.75, 0.81, r'$y=x$', rotation=38, ha='center', va='center',
    #          color='darkslategray', alpha=0.65)
    # plt.text(tx1, ty1, 'A', transform=plt.gcf().transFigure)
    if loo:
        f_save = figdir + '/fig5.{}.leave-one-out.png'.format(fldr)
    else:
        f_save = figdir + '/fig5.{}.png'.format(fldr)

    plt.savefig(f_save, dpi=512)
    plt.close()


# fldr = 'cadd94_gmask_mnb_378'
# fldr = 'cadd94_gmask_v1.6_without_bstat'
# fldr = 'mockYRI'
# fldr = 'dropped_chrom_results'
fldr = 'cadd94_gmask_v1.6_without_bstat'
# fldr = 'cadd94_gmask_v1.6'

figure_5_final_format(fldr, loo=True)
figure_5_final_format(fldr)


#%% COLLATED PLOT AND MULTIPLE SORTED PREDICTIONS COMBINED
def collate_and_sort(flist, llist, clist, lwdths, sname, fletters=('d', 'e'),
                     legend_on=False):

    obs_flag = False
    # plt.figure(figsize=(6.5, 2))
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.062, bottom=0.152, right=0.99,
                        top=0.98, wspace=0.25)

    # 1. COLLATED PLOT
    # plt.subplot(121)
    ax1 = plt.subplot(1, 3, (1,2))
    # format_panels(ax1)
    for (i, rdir) in enumerate(flist):
        f_ape = root_dir + '/result/final_files/{}/nonsyn.collate.5.00e-03width.npy'.format(
            rdir)
        f_mcv = root_dir + '/result/collate/YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy'
        bins, div, pi, pr, cnts = np.load(f_ape).T
        # _, _, _, mcv, mcv_cnts = np.load(f_mcv).T
        # bins, pi, div, pred, cnts = np.load(f_in).T
        obs = pi / cnts
        prd = pr / cnts
        # mcv /= mcv_cnts
        # y_max = pi0
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        # mcv /= y_max
        # plt.title(rdir)
        # pc = '{}%'.format(100-int(rdir.split('_')[1][4:]))
        if not obs_flag:
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=2)
            obs_flag = True
        # plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=1.1)
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.8, lw=lwdths[i])

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=1)
    plt.xticks(y=0.04)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=1)
    plt.yticks(x=0.02)
    if legend_on:
        plt.legend(loc='lower right', ncol=1, frameon=1,
                   framealpha=0.75, facecolor='white', handlelength=0.8,
                   borderaxespad=0.3, columnspacing=0.8, prop=dict(size=9))
    # plot inset
    obs_flag = False
    axins = inset_axes(ax1, width="100%", height="100%", loc=6,
                       bbox_to_anchor=(0.021, 0.1, 0.35, 0.35),
                       bbox_transform=ax1.transAxes)
    format_panels(axins)
    for (i, fl) in enumerate(flist):
        f_ape = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(fl)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        msk = (bins >= -0.05) & (bins <= 0.05)
        obs = pi / cnts
        prd = pr / cnts
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        if not obs_flag:
            plt.plot(bins[msk], obs[msk], label='observed',
                     color='darkslategray', lw=1.5)
            obs_flag = True
        plt.plot(bins[msk], prd[msk], label=llist[i], color=clist[i], alpha=0.9,
                 lw=lwdths[i], ls='-')
    plt.xticks([-0.05, 0.05], y=0.07, fontsize=8)
    plt.xlim(-0.05, 0.05)
    axins.set_yticks([])
    plt.ylim(0.81, 1.05)

    # 2. PREDICTED VS. OBSERVED PLOT
    # plt.subplot(122)
    ymin, xmin = 1, 1
    xmax, ymax = 1, 1
    ax2 = plt.subplot(133)
    # format_panels(ax2)
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='k',
             ls='--', alpha=1)
    for (i, rdir) in enumerate(flist):
        span = 0.1
        fldr = flist[i]
        color = clist[i]
        # get 100 points data
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

        # plot horizontal line at y=1
        plt.axhline(y=1, color='k', alpha=0.8, ls='-')

        # plot predicted vs. observed
        plt.plot(pred, pi, marker='o', ms=2, color='white', lw=0)
        plt.plot(pred, pi, marker='o', ms=2, color=color, lw=0, alpha=0.5)

        # # plot LOESS line
        # # plt.plot(prlo, pilo, lw=2, color='white')
        # plt.plot(prlo, pilo, lw=lwdths[i], color=color, alpha=0.8)
        # xmin = min(xmin, prlo.min(), pred.min())
        # xmax = max(xmax, prlo.max(), pred.max())
        # ymin = min(ymin, pilo.min(), pi.min())
        # ymax = max(ymax, pilo.max(), pi.max())

        for (i, rdir) in enumerate(flist):
            span = 0.1
            fldr = flist[i]
            color = clist[i]
            # get 100 points data
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
            # plot LOESS line
            # plt.plot(prlo, pilo, lw=2, color='white')
            plt.plot(prlo, pilo, lw=lwdths[i], color=color, alpha=0.8)
            xmin = min(xmin, prlo.min(), pred.min())
            xmax = max(xmax, prlo.max(), pred.max())
            ymin = min(ymin, pilo.min(), pi.min())
            ymax = max(ymax, pilo.max(), pi.max())


        # if color == 'darkorange':
        #     plt.plot(prlo, pilo, lw=1.5, color=color)
    print('ylims {} {}'.format(ymin, ymax))
    print('xlims {} {}'.format(xmin, xmax))
    if 'McVicker' in llist:
        plt.ylim(0., 1.15)
        plt.xlim(0., 1.02)
        xtick = np.arange(0.1, 1.01, 0.1)
        ytick = np.arange(0.1, 1.11, 0.1)
    else:
        axmin = 0.48
        tckmin = 0.5
        # plt.ylim(0.48, 1.15)
        # plt.xlim(0.48, 1.02)
        # xtick = np.arange(0.5, 1.01, 0.1)
        # ytick = np.arange(0.5, 1.11, 0.1)
        plt.ylim(ymin-0.02, ymax+0.02)
        plt.xlim(xmin-0.02, xmax+0.02)
        xtick = np.arange(0.1, 1.01, 0.1)
        ytick = np.arange(0.1, 1.21, 0.1)

    plt.xticks(xtick, y=0.04)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=1)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=1)
    plt.yticks(ytick, x=0.04)
    plt.ylim(ymin - 0.02, ymax + 0.02)
    plt.xlim(xmin - 0.02, xmax + 0.02)
    plt.text(xmin, 1.03, 'without linked selection', ha='left',
             va='center', fontsize=10)
    # solve y=x rotation
    adlen = xmax+0.02-(ymin-0.02)
    oplen = adlen * (xmax-xmin+0.04) / (ymax-ymin+0.04)
    rot = np.arctan((oplen/adlen)) * (180.0 / np.pi)
    # rot = 45 * oplen/adlen
    print('rotation = {}'.format(rot))
    plt.text(0.76, 0.64, r'$y=x$', rotation=rot, ha='left', va='bottom',
             color='k', alpha=1)
    plt.text(0.065, 0.915, fletters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.73, 0.915, fletters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/mainfigs/{}.collated.sort.combined.png'.format(sname)
    plt.savefig(f_save, dpi=512)
    plt.close()

flist = ['cadd94_gmask_v1.6_without_bstat', 'cadd94_gmask_v1.6_without_bstat_jackknife_results',
         'dropped_chrom_results']
llist = ['all data', 'leave 2 Mb out', 'leave autosome out']
clist = ['darkorange', 'deepskyblue', 'fuchsia']
sname = 'compare_LOO'

# rsq_3_size(flist, llist2, sname, xlab, rotation=0)
# combined_chr1(flist, llist, clist, sname)
collate_and_sort(flist, llist, clist, lwidths, sname, fletters=('c', 'd'))


#%%