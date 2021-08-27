__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, cst_from_fldr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from figures.common_functions import format_panels, predict_loess


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


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


#%% COLLATED FOR DIFFERENT POPS
def collated_ns_other_pops(fldr, fcl, xlab, color, title, letter):
    # fixed variables:
    # color = 'darkorange'
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
    plt.title(title, y=0.85)
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
        plt.ylim(0.8, 1.165)
        plt.yticks(x=0.03)
    plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
               borderpad=0.2, framealpha=0.75, facecolor='white',
               handlelength=1.2, labelspacing=0.25)
    if letter is not None:
        plt.text(0.145, 0.915, letter, transform=plt.gcf().transFigure,
                 fontweight='bold')

    f_save = final_dir + '/sfigs/{}.{}.collated.png'.format(fldr, fcl)
    plt.savefig(f_save, dpi=512)
    plt.close()


fcl = 'nonsyn'
xlab = 'distance to nearest NS substitution (cM)'
letters = 'a b c d'.split()
titles = ['European data from CEU', 'South-Asian data from GIH',
          'East-Asian data from JPT', 'Amerindian data from MXL']
use_pop = ['CEU', 'GIH', 'JPT', 'MXL']
cols = ['steelblue', 'darkturquoise', 'purple', 'fuchsia']
for i in xrange(4):
    pop = use_pop[i]
    fldr = '{}_cadd94_gmask_v1_6_without_bstat'.format(pop)
    collated_ns_other_pops(fldr, fcl, xlab, cols[i], titles[i], letters[i])


#%% COLLATED PLOTS AROUND OTHER FOCAL SITES
def collate_other_focals(fldr, fcl, xlab, letter=None):
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

    plt.figure(figsize=(3.25, 2.5))
    plt.subplots_adjust(left=0.135, bottom=0.2, right=0.995, top=0.995)
    ax1 = plt.subplot(111)
    format_panels(ax1)

    plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(bins, prd, label=label, color=color, alpha=0.8, lw=1.1)
    plt.xlabel(xlab, labelpad=3)
    plt.xticks(y=0.04)
    # plt.xlim(-0.41, 0.41)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=1)
    if fcl == 'YRI_syn_s1':
        plt.ylim(0.81, 1.13)
    if fcl == 'hc_derived_cons':
        plt.yticks([0.7, 1.0, 1.3], x=0.03)
    else:
        plt.yticks(x=0.03)
    if letter == 'a':
        plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
                   borderpad=0.2, framealpha=0.75, facecolor='white',
                   handlelength=1.2, labelspacing=0.25)
    if letter is not None:
        plt.text(0.145, 0.935, letter, transform=plt.gcf().transFigure,
                 fontweight='bold')

    f_save = final_dir + '/sfigs/{}.{}.collated.png'.format(fldr, fcl)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat'
focals = ['YRI_syn_s1', 'hc_derived_cons', 'exon', 'fish_cons94_new_exonic']
xlabels = ['synonymous\nsubstitution', 'conserved\nsubstitution', 'exon',
           'conserved\nexonic segment']
letters = 'a b c d'.split()
idx = 0
indices = range(4)
for idx in indices:
    xlab = 'distance to nearest {} (cM)'.format(xlabels[idx])
    fcl = focals[idx]
    collate_other_focals(fldr, fcl, xlab, letters[idx])


#%% OBS/PRED FOR DIFFERENT POPS
def figure_5_inset(fldr, pct, span, color, lcolor, letter=None):
    # fixed variables
    num = 100
    axmin, axmax = 0.55, 1.15

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
    idx = np.where(ipred >= pct)[0][0]
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
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    # BOX ON INSET DATA
    inset_lw = 0.5
    rect = Rectangle(xy=(pct, pct), width=0.1, height=0.24,
                     linewidth=inset_lw, edgecolor='red',
                     facecolor='none')
    ax1.add_patch(rect)

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
    axins.plot(ipred, pr_loess, label=loess_lab, color=lcolor,
               alpha=0.95)
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
    if letter is not None:
        plt.text(0.16, 0.935, letter, transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/{}.sort.inset.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


letters = 'A B C D E'.split()
use_pop = ['CEU', 'GIH', 'JPT', 'MXL']
cols = ['steelblue', 'darkturquoise', 'purple', 'fuchsia']
lcols = ['mediumvioletred', 'mediumvioletred', 'turquoise', 'turquoise']
for i in [0, 1]:
    pop = use_pop[i]
    fldr = '{}_cadd94_gmask'.format(pop)
    figure_5_inset(fldr, 0.9, 0.1, cols[i], lcols[i], letters[i])


#%% OBSERVED VS. PREDICTED FOR DIFFERENT POPS (USED IN PAPER)
def figure_5_diff_pops(fldr, color, lcolor, title, letter, span=0.1):
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

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)
    plt.title(title, y=0.91)

    axmin, axmax = 0.53, 1.25
    tx1, tx2 = 0.16, 0.61
    ty1, ty2 = 0.94, 0.48
    xtick = np.arange(0.5, 1.2, 0.1)

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='k', ls='--')
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    # plot predicted vs. observed
    plt.plot(pred, pi, marker='o', ms=5, color=color, lw=0,
             alpha=0.5)
    # plot LOESS line
    plt.plot(prlo, pilo, lw=2, color=lcolor)
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.21, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xticks(xtick, y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.45, axmax)
    plt.xlim(axmin, 1.02)
    # solve y=x rotation
    # adlen = xmax+0.02-(ymin-0.02)
    # oplen = adlen * (xmax-xmin+0.04) / (ymax-ymin+0.04)
    adlen = 1.02 - axmin
    oplen = (1.02 - axmin) * ((1.02-axmin) / (axmax-0.45))
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')
    plt.text(tx1, ty1, letter, transform=plt.gcf().transFigure,
             fontweight='bold')
    f_save = final_dir + '/sfigs/{}.sort.inset.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


letters = 'a b c d'.split()
titles = ['European data from CEU', 'South-Asian data from GIH',
          'East-Asian data from JPT', 'Amerindian data from MXL']
use_pop = ['CEU', 'GIH', 'JPT', 'MXL']
cols = ['steelblue', 'darkturquoise', 'purple', 'fuchsia']
lcols = ['mediumblue', 'DarkCyan', 'indigo', 'MediumPurple']
for i in [0, 1, 2, 3]:
    pop = use_pop[i]
    fldr = '{}_cadd94_gmask_v1_6_without_bstat'.format(pop)
    figure_5_diff_pops(fldr, cols[i], lcols[i], titles[i], letters[i])


#%% CHR1 RESULTS ACROSS POPULATIONS (USED IN PAPER)
def combined_chr1(fldr, color, title, fletter):
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/sfigs'
    obs_flag = False
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.07, bottom=0.155, right=0.995, top=0.995)
    xi = None
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title(title, y=0.85)
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    chlist = [f for f in os.listdir(fdir) if ('chr1' in f) and
             (f.endswith('.txt')) and 'TMRCA' not in f]
    f_name = fdir + chlist[0]
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    obs[(obs < 0.1) | (obs > 2)] = np.nan

    plt.plot(xi, obs, label='observed', color='darkslategray', lw=2)
    plt.plot(xi, prd, color=color, lw=1.5, alpha=0.8, label='predicted')

    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5, 2], x=0.0125)
    # plt.ylim(-0.4, 1.9)
    plt.ylim(0, 2.2)
    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x / 25) % 2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])
    loc = 'lower right'
    plt.legend(loc=loc, ncol=2, frameon=1,
               framealpha=0.75, facecolor='white', handlelength=0.8,
               borderaxespad=0.3, columnspacing=0.8)

    # figure letter
    plt.text(0.075, 0.92, fletter, transform=plt.gcf().transFigure,
             fontweight='bold')
    f_save = sdir + '/{}.chr1.combined.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


letters = 'a b c d'.split()
titles = ['European data from CEU', 'South-Asian data from GIH',
          'East-Asian data from JPT', 'Amerindian data from MXL']
use_pop = ['CEU', 'GIH', 'JPT', 'MXL']
cols = ['steelblue', 'darkturquoise', 'purple', 'fuchsia']
for i in xrange(4):
    pop = use_pop[i]
    fldr = '{}_cadd94_gmask_v1_6_without_bstat'.format(pop)
    combined_chr1(fldr, cols[i], titles[i], letters[i])


#%% COLLATED PLOTS FOR REP. POPS. USING YRI MAP AND SELF-MAPS
def collated_plot_compare_to_yri_pred(fldr, ttl, ltr, clr, sname):
    obs_flag = False
    # llist = ['observed', 'predicted', 'YRI-map predicted']
    plt.figure(figsize=(4, 2.5))
    plt.subplots_adjust(left=0.145, bottom=0.15, right=1, top=1)

    # COLLATED PLOT
    ax1 = plt.subplot(111)
    format_panels(ax1)
    # plt.text(0, 1.08, ttl, ha='center', va='center', color='k')
    plt.title(ttl, y=0.85)
    fdir = final_dir + '/{}/'.format(fldr)
    f_ape = fdir + 'nonsyn.collate.5.00e-03width.npy'
    bins, div, pi, pr, cnts = np.load(f_ape).T
    obs = pi / cnts
    prd = pr / cnts
    y_max = obs.mean()
    obs /= y_max
    prd /= y_max
    bins += 2.5e-3
    # plot observed diversity
    plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.5)
    # plot predicted diversity based on diversity data used for inference
    plt.plot(bins, prd, label='predicted', color=clr, alpha=0.7, lw=0.8,
             ls='-')
    # plot predicted diversity based on YRI diversity data (rescaled)
    f_yri = fdir + 'YRI_sorted.nonsyn.collate.5.00e-03width.npy'
    bins, div, pi, pr, cnts = np.load(f_yri).T
    obs = pi / cnts
    prd = pr / cnts
    y_max = obs.mean()
    obs /= y_max
    prd /= y_max
    bins += 2.5e-3
    plt.plot(bins, prd, label='YRI-predicted', color=clr, alpha=0.7, lw=0.8,
             ls='--')

    plt.xlabel('distance to nearest amino acid substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    plt.ylim(0.81, 1.12)
    plt.legend(prop=dict(size=7), loc='lower right', handletextpad=0.5,
               borderaxespad=0.05, frameon=True, framealpha=0.7,
               facecolor='white')

    # zoom on middle
    obs_flag = False
    axins = inset_axes(ax1, width="100%", height="100%", loc=6,
                       bbox_to_anchor=(0.02, 0.1, 0.35, 0.35),
                       bbox_transform=ax1.transAxes)
    format_panels(axins)
    sort_type = ['', 'YRI_sorted.']
    line_type = ['-', '--']
    lbls = ['predicted', 'YRI-map predicted']
    for (i, st) in enumerate(sort_type):
        f_ape = fdir + st + 'nonsyn.collate.5.00e-03width.npy'
        bins, div, pi, pr, cnts = np.load(f_ape).T
        bins += 2.5e-3
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
        plt.plot(bins[msk], prd[msk], label=lbls[i], color=clr, alpha=0.7,
                 lw=0.8, ls=line_type[i])

    plt.xticks([-0.05, 0.05], y=0.07, fontsize=8)
    plt.xlim(-0.05, 0.05)
    axins.set_yticks([])
    plt.ylim(0.81, 1.05)

    f_save = final_dir + '/sfigs/{}.collated.png'.format(sname)
    plt.savefig(f_save, dpi=512)
    plt.close()


folder = 'JPT_cadd94_gmask'
color = 'purple'
sname = 'JPT_collate_self_and_YRI_sort'
colors = ['purple', 'darkturquoise', 'fuchsia', 'steelblue']
pops = 'JPT GIH MXL CEU'.split()
for i in xrange(4):
    folder = pops[i] + '_cadd94_gmask'
    sname = pops[i] + '_collate_self_and_YRI_sort'
    collated_plot_compare_to_yri_pred(folder, pops[i], colors[i], sname)


#%% COLLATED PLOTS FOR REP. POPS. USING YRI MAP AND SELF-MAPS -- data only
def collated_plot_dataonly(flist, llist, clist, sname):
    plt.figure(figsize=(4, 2.5))
    plt.subplots_adjust(left=0.145, bottom=0.15, right=1, top=1)

    # COLLATED PLOT
    ax1 = plt.subplot(111)
    format_panels(ax1)
    for (i, fl) in enumerate(flist):
        f_ape = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(fl)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        bins += 2.5e-3
        obs = pi / cnts
        y_max = obs.mean()
        obs /= y_max
        plt.plot(bins, obs, label=llist[i], color=clist[i], lw=0.75, alpha=0.75)

    plt.xlabel('distance to nearest amino acid substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    plt.ylim(0.81, 1.12)
    plt.legend(prop=dict(size=7), loc='lower right', handletextpad=0.5,
               borderaxespad=0.05, frameon=True, framealpha=0.7,
               facecolor='white')

    # zoom on middle
    axins = inset_axes(ax1, width="100%", height="100%", loc=6,
                       bbox_to_anchor=(0.02, 0.1, 0.35, 0.35),
                       bbox_transform=ax1.transAxes)
    format_panels(axins)
    for (i, fl) in enumerate(flist):
        f_ape = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(fl)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        bins += 2.5e-3
        msk = (bins >= -0.05) & (bins <= 0.05)
        obs = pi / cnts
        prd = pr / cnts
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        plt.plot(bins[msk], obs[msk], label=llist[i], color=clist[i], alpha=0.7,
                 lw=0.8, ls='-')

    plt.xticks([-0.05, 0.05], y=0.07, fontsize=8)
    plt.xlim(-0.05, 0.05)
    axins.set_yticks([])
    plt.ylim(0.81, 1.05)

    f_save = final_dir + '/sfigs/{}.collated.png'.format(sname)
    plt.savefig(f_save, dpi=512)
    plt.close()

colors = 'darkorange steelblue darkturquoise purple fuchsia'.split()
pops = 'YRI CEU GIH JPT MXL'.split()
folders = ['cadd94_gmask_mnb_378'] + ['{}_cadd94_gmask'.format(p) for p in pops[1:]]

collated_plot_dataonly(folders, pops, colors, 'continental_pops_dataonly')

#%%