__author__ = 'davidmurphy'

import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir, \
    human_autosomes, RunStruct, cst_from_fldr
from figures.common_functions import format_panels, final_dir, \
    rsq_three_scales, get_sortplot_rsq, predict_loess
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr


llh_thresholds = [0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8]
pre_thresholds = [t**(7.4/0.995) for t in llh_thresholds]


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


#%% GATHER RSQ AND OTHER DATA FROM FINAL DIRS
cad_bth_rsq = []
cad_mnb_rsq = []

# lists for sort rsq values
# ape_bth_srsq = []
# ape_mnb_srsq = []
cad_bth_srsq = []
cad_mnb_srsq = []

# lists for utot values
# ape_bth_utot = []
# ape_mnb_utot = []
cad_bth_utot = []
cad_mnb_utot = []

# folder formats the different thresholds
# ape_bth_fmt = 'fish_cons94_gmask_bth_{}'
# ape_mnb_fmt = 'fish_cons94_gmask_mnb_{}'
# ape_bth_fmt = 'ape_cons94_gmask_bth_{}'
# ape_mnb_fmt = 'ape_cons94_gmask_mnb_{}'
cad_bth_fmt = 'cadd94_gmask_v1.6_without_bstat_bth_{}'
cad_mnb_fmt = 'cadd94_gmask_v1.6_without_bstat_minb_{}'

cadd_b0 = rsq_three_scales(cad_bth_fmt.format('000'))
cad_bth_rsq.append(cadd_b0)
cad_mnb_rsq.append(cadd_b0)

# ape_b0 = rsq_three_scales(ape_bth_fmt.format('000'))
# ape_bth_rsq.append(cadd_b0)
# ape_mnb_rsq.append(cadd_b0)

# ape_bth_utot.append(cst_from_fldr(ape_bth_fmt.format('000')).stat.utot[0])
# ape_mnb_utot.append(cst_from_fldr(ape_bth_fmt.format('000')).stat.utot[0])
# cad_bth_utot.append(cst_from_fldr(cad_bth_fmt.format('000')).stat.utot[0])
# cad_mnb_utot.append(cst_from_fldr(cad_bth_fmt.format('000')).stat.utot[0])
#
# ape_bth_srsq.append(get_sortplot_rsq(ape_bth_fmt.format('000'), 100))
# ape_mnb_srsq.append(get_sortplot_rsq(ape_bth_fmt.format('000'), 100))
# cad_bth_srsq.append(get_sortplot_rsq(cad_bth_fmt.format('000'), 100))
# cad_mnb_srsq.append(get_sortplot_rsq(cad_bth_fmt.format('000'), 100))
# 105, 121, 140
# the strings used in each threshold folder
# 140 121 105 913 797 695 603 520 444 310 250 194
# mnb_tkn = '105 121 140 194 250 310 375 444 520 603 695 797 913'.split()
# mnb_tkn = '140 121 105 913 797 695 603 520 444 375 310 250 194'.split()
# bth_tkn = '200 250 300 350 400 450 500 550 600 650 700 750 800'.split()

# mnb_tkn = '140 121 105 913 797 695 603 520 444 375 310 194'.split()
# bth_tkn = '200 250 300 350 400 450 500 550 600 650 700 800'.split()
mnb_all = '119 891 678 591 513 442 378 319 264 165'.split()
bth_all = '200 300 400 450 500 550 600 650 700 800'.split()
# get data from each folder
for i in xrange(len(mnb_all)):
    # # get folders for each type of threshold
    # ape_bth_fldr = ape_bth_fmt.format(bth_all[i])
    # ape_mnb_fldr = ape_mnb_fmt.format(mnb_all[i])
    cad_bth_fldr = cad_bth_fmt.format(bth_all[i])
    cad_mnb_fldr = cad_mnb_fmt.format(mnb_all[i])

    # # replace min b of -4.44 with standard -4.61
    # if mnb_all[i] == '444':
    #     ape_mnb_fldr = 'ape94_gmask_01'
    #     cad_mnb_fldr = 'cadd93_gmask_01'

    # get all of the R^2 values
    # ape_bth_rsq.append(rsq_three_scales(ape_bth_fldr))
    # ape_mnb_rsq.append(rsq_three_scales(ape_mnb_fldr))
    cad_bth_rsq.append(rsq_three_scales(cad_bth_fldr))
    cad_mnb_rsq.append(rsq_three_scales(cad_mnb_fldr))

    # # get all of the sort map R^2 values
    # ape_bth_srsq.append(get_sortplot_rsq(ape_bth_fldr, 100))
    # ape_mnb_srsq.append(get_sortplot_rsq(ape_mnb_fldr, 100))
    # cad_bth_srsq.append(get_sortplot_rsq(cad_bth_fldr, 100))
    # cad_mnb_srsq.append(get_sortplot_rsq(cad_mnb_fldr, 100))
    #
    # # get all of the utot values
    # ape_bth_utot.append(cst_from_fldr(ape_bth_fldr).stat.utot[0])
    # ape_mnb_utot.append(cst_from_fldr(ape_mnb_fldr).stat.utot[0])
    # cad_bth_utot.append(cst_from_fldr(cad_bth_fldr).stat.utot[0])
    # cad_mnb_utot.append(cst_from_fldr(cad_mnb_fldr).stat.utot[0])

# convert all lists to arrays
# ape_bth_rsq = np.array(ape_bth_rsq)
# ape_mnb_rsq = np.array(ape_mnb_rsq)

# ape_bth_srsq = np.array(ape_bth_srsq)
# ape_mnb_srsq = np.array(ape_mnb_srsq)
# cad_bth_srsq = np.array(cad_bth_srsq)
# cad_mnb_srsq = np.array(cad_mnb_srsq)
#
# ape_bth_utot = np.array(ape_bth_utot)
# ape_mnb_utot = np.array(ape_mnb_utot)
# cad_bth_utot = np.array(cad_bth_utot)
# cad_mnb_utot = np.array(cad_mnb_utot)
cad_bth_rsq = np.array(cad_bth_rsq)
cad_mnb_rsq = np.array(cad_mnb_rsq)


#%% R^2 OVER 3 ORDERS OF MAGNITUDE
def threshold_rsq_plot(anno):
    # format fig
    plt.figure(figsize=(6.5, 2.16))

    # bth_tkn = '000 200 250 300 350 400 450 500 550 600 650 700 800'.split()
    bth_tkn = bth_all

    xi = np.arange(len(bth_tkn)+1)
    xtcks = [0]+['{:.2f}'.format(0.001 * float(b)) for b in bth_tkn]

    # PANEL 1: 10KB
    ax1 = plt.subplot(131)
    format_panels(ax1)
    plt.subplots_adjust(bottom=0.22, left=0.1, right=0.995, top=0.9,
                        wspace=0.25)
    yttl = 0.96

    plt.plot(xi, cad_mnb_rsq[:, 0], marker='o', label='lookup',
             color='dodgerblue', alpha=0.7, ms=3, lw=0.8)
    plt.plot(xi, cad_bth_rsq[:, 0], marker='s', label='maximization',
             color='forestgreen', alpha=0.7, ms=3, lw=0.8)
    plt.xticks(xi, xtcks, rotation=90, y=0.03)
    plt.xlabel(r'$B$ threshold', labelpad=2)
    ytck = np.arange(0.135, 0.16, 0.005)
    plt.yticks(ytck, x=0.05)
    plt.ylabel(r'variance explained ($R^2$)', labelpad=1)
    plt.title('10 kb windows', y=yttl)


    plt.legend(frameon=True, framealpha=0.75, facecolor='white',
               handlelength=0.8, borderaxespad=0.3, loc='lower center')

    # PANEL 2: 100KB
    ax2 = plt.subplot(132)
    format_panels(ax2)
    # plt.plot(xi, ape_mnb_rsq[:,1], marker='o', label='ape 100Kb precalc',
    #          color='darkorange', alpha=0.7, ms=3, lw=0.8)
    # plt.plot(xi, ape_bth_rsq[:,1], marker='s', label='ape 100Kb likelihood',
    #          color='firebrick', alpha=0.7, ms=3, lw=0.8)
    plt.plot(xi, cad_mnb_rsq[:,1], marker='o', label='CADD (pc)',
             color='dodgerblue', alpha=0.7, ms=3, lw=0.8)
    plt.plot(xi, cad_bth_rsq[:,1], marker='s', label='CADD 100Kb likelihood',
             color='forestgreen', alpha=0.7, ms=3, lw=0.8)
    plt.xticks(xi, xtcks, rotation=90, y=0.03)
    plt.xlabel(r'$B$ threshold', labelpad=2)
    ytck = np.arange(0.35, 0.41, 0.01)
    plt.yticks(ytck, x=0.05)
    # plt.ylabel(r'$\mathrm{R^2}$ 100kb', labelpad=1)
    plt.title('100 kb windows', y=yttl)


    # PANEL 3: 1MB
    ax3 = plt.subplot(133)
    format_panels(ax3)
    # plt.plot(xi, ape_mnb_rsq[:, 2], marker='o', label='ape 1Mb precalc',
    #          color='darkorange', alpha=0.7, ms=3, lw=0.8)
    # plt.plot(xi, ape_bth_rsq[:, 2], marker='s', label='ape 1Mb likelihood',
    #          color='firebrick', alpha=0.7, ms=3, lw=0.8)
    plt.plot(xi, cad_mnb_rsq[:, 2], marker='o', label='CADD 1Mb precalc',
             color='dodgerblue', alpha=0.7, ms=3, lw=0.8)
    plt.plot(xi, cad_bth_rsq[:, 2], marker='s', label='CADD 1Mb likelihood',
             color='forestgreen', alpha=0.7, ms=3, lw=0.8)
    plt.xticks(xi, xtcks, rotation=90, y=0.03)
    plt.xlabel(r'$B$ threshold', labelpad=2)
    ytck = np.arange(0.54, 0.61, 0.01)
    plt.yticks(ytck, x=0.05)
    # plt.ylabel(r'$\mathrm{R^2}$ 1Mb', labelpad=1)
    plt.title('1 Mb windows', y=yttl)

    f_save = final_dir + '/sfigs/fig_S8_rsq_compare_CADD_{}_B_thresholds.png'.format(anno)
    plt.savefig(f_save, dpi=512)
    plt.close()


threshold_rsq_plot('cadd_alone')



#%% OBSERVED VS. PREDICTED ACROSS THRESHOLDS (UPDATED JUNE 2021)
def figure_5_diff_thresholds(flist, llist, clist, lclist, span=0.1):
    axmin, axmax = 0.32, 1.4
    tx1, tx2 = 0.16, 0.61
    ty1, ty2 = 0.945, 0.48
    xtick = np.arange(0, 1.4, 0.1)

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)
    # plot y=x line label=r'$y=x$',
    plt.plot([axmin, axmax], [axmin, axmax],
             color='k', ls='--')

    # plot each data set
    for i in xrange(len(flist)):
        fldr = flist[i]
        color = clist[i]
        lcolor = lclist[i]
        lbl = llist[i]
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

        # plot LOESS line
        # plt.plot(prlo, pilo, lw=2, color='white')
        plt.plot(prlo, pilo, lw=1.5, color=color, alpha=0.8, label=lbl)
        if color == 'darkorange':
            plt.plot(prlo, pilo, lw=1.5, color=color)
    plt.legend(loc='upper center', ncol=2, handletextpad=0.5,
               borderaxespad=0.2, frameon=True, framealpha=0.9,
               facecolor='white', title='Threshold B-value')
    # plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
    #          va='center', fontsize=11)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.4, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xticks(xtick, y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.45, axmax)
    plt.xlim(axmin, 1.02)
    # plt.ylim(0.55, 1.2)
    # plt.xlim(0.55, 1.02)
    # solve y=x rotation
    adlen = axmax - 0.52
    oplen = adlen * (1.02 - axmin) / (axmax - 0.52)
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')

    plt.text(tx1, ty1, 'a', transform=plt.gcf().transFigure, fontweight='bold')
    f_save = final_dir + '/sfigs/fig_S7a.obs.pred.across.threshold.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


# mnb_vals = '119 512 378'.split()
bth_vals = [0.2, 0.5, 0.6]
mnblist = ['cadd94_gmask_v1.6_without_bstat_bth_000',
           'cadd94_gmask_v1.6_without_bstat_minb_119',
           'cadd94_gmask_v1.6_without_bstat_minb_513',
           'cadd94_gmask_v1.6_without_bstat']
llist = ['none'] + [r'B = {:.2f}'.format(b) for b in 0.2, 0.5, 0.6]
clist = 'Maroon dodgerblue deeppink darkorange'.split()
# llist = ['old filtering', 'new filtering']
# twofish = ['fish_cons94_gmask_mnb_378', 'fish_cons94_new']
figure_5_diff_thresholds(mnblist, llist, clist, clist)


#%% COLLATED PLOT RESULTS (COMBINED) - UPDATED JUNE 2021
def collated_plot(flist, llist, clist):
    obs_flag = False
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(left=0.16, bottom=0.125, right=0.995, top=0.99)

    # COLLATED PLOT
    ax1 = plt.subplot(111)
    format_panels(ax1)
    for (i, fl) in enumerate(flist):
        f_ape = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(fl)
        bins, div, pi, pr, cnts = np.load(f_ape).T
        obs = pi / cnts
        prd = pr / cnts
        y_max = obs.mean()
        obs /= y_max
        prd /= y_max
        if not obs_flag:
            plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.5)
            obs_flag = True
        plt.plot(bins, prd, label=llist[i], color=clist[i], alpha=0.7, lw=0.8,
                 ls='-')

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=3)
    plt.xticks(y=0.02)
    plt.xlim(-0.22, 0.22)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=3)
    plt.yticks(x=0.02)
    plt.ylim(0.82, 1.0525)
    plt.legend(prop=dict(size=10), loc='lower right', handletextpad=0.5,
               borderaxespad=0.2, frameon=True, framealpha=0.9,
               facecolor='white')
    plt.text(0.169, 0.949, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    # zoom on middle
    obs_flag = False
    # bbox_to_anchor = x0, y0, x1, y1
    axins = inset_axes(ax1, width="100%", height="100%", loc=6,
                       bbox_to_anchor=(0.03, 0.06, 0.35, 0.35),
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
        plt.plot(bins[msk], prd[msk], label=llist[i], color=clist[i], alpha=0.7,
                 lw=0.8, ls='-')

    plt.xticks([-0.05, 0.05], y=0.07, fontsize=8)
    plt.xlim(-0.05, 0.05)
    axins.set_yticks([])
    plt.ylim(0.81, 1.049)

    f_save = final_dir + '/sfigs/fig_S7b.collated.across.threshold.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


mnblist = ['cadd94_gmask_v1.6_without_bstat_bth_000',
           'cadd94_gmask_v1.6_without_bstat_minb_119',
           'cadd94_gmask_v1.6_without_bstat_minb_513',
           'cadd94_gmask_v1.6_without_bstat']
llist = ['none'] + [r'B = {:.2f}'.format(b) for b in 0.2, 0.5, 0.6]
clist = 'Maroon dodgerblue deeppink darkorange'.split()
collated_plot(mnblist, llist, clist)


#%% COMRPARE MULTIPLE CHR1 PREDICTIONS
def combined_chr1(flist, llist, clist):
    # make a folder using sname to save files in
    sdir = root_dir + '/result/final_files/sfigs'
    obs_flag = False
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.07, bottom=0.155, right=0.995, top=0.995)
    xi = None
    ax = plt.subplot(111)
    format_panels(ax)
    for (i, rdir) in enumerate(flist):
        fdir = root_dir + '/result/final_files/{}/'.format(rdir)
        chlist = [f for f in os.listdir(fdir) if ('chr1' in f) and
                 (f.endswith('.txt')) and 'TMRCA' not in f]
        f_name = fdir + chlist[0]
        prd, obs, num = np.loadtxt(f_name).T
        xi = np.arange(0, prd.size / 2.0, 0.5)
        pi_mean = np.nanmean(obs)
        prd /= pi_mean
        obs /= pi_mean
        obs[(obs < 0.25) | (obs > 1.75)] = np.nan

        if not obs_flag:
            plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
            obs_flag = True

        plt.plot(xi, prd, color=clist[i], lw=0.8, alpha=0.8, label=llist[i])

    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(-0.4, 1.9)
    plt.ylim(0, 1.9)
    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x / 25) % 2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])
    plt.legend(loc='lower center', ncol=len(llist)+1, frameon=1,
               framealpha=0.75, facecolor='white', handlelength=0.8,
               borderaxespad=0.3, columnspacing=0.8)
    inset_lw = 0.5
    rect = Rectangle(xy=(45, 0.31), width=10, height=0.8,
                     linewidth=inset_lw, edgecolor='red',
                     facecolor='none')
    ax.add_patch(rect)
    # figure letter
    # plt.text(0.01, 0.93, fletter, transform=plt.gcf().transFigure)
    f_save = sdir + '/fig_S9.chr1.across.thresholds.png'
    plt.savefig(f_save, dpi=256)
    plt.close()


mnblist = ['cadd94_gmask_v1.6_without_bstat_bth_000',
           'cadd94_gmask_v1.6_without_bstat_minb_119',
           'cadd94_gmask_v1.6_without_bstat_minb_513',
           'cadd94_gmask_v1.6_without_bstat']
llist = ['none'] + [r'B = {:.2f}'.format(b) for b in 0.2, 0.5, 0.6]
clist = 'Maroon dodgerblue deeppink darkorange'.split()
combined_chr1(mnblist, llist, clist)


#%% CHECK PEARSON ON THE MIDDLE 90% OF VALUES FOR THRESHOLDS
def check_calibration_pearson(flist, llist):
    for i in xrange(len(flist)):
        fldr = flist[i]
        labl = llist[i]
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 100)
        div, pi, pred = np.loadtxt(sort_file).T
        pcorr = pearsonr(pi[5:95], pred[5:95])
        print labl, pcorr


mnblist = ['cadd94_gmask_v1.6_without_bstat_bth_000',
           'cadd94_gmask_v1.6_without_bstat_minb_119',
           'cadd94_gmask_v1.6_without_bstat_minb_513',
           'cadd94_gmask_v1.6_without_bstat']
llist = ['none'] + [r'B = {:.2f}'.format(b) for b in 0.2, 0.5, 0.6]
check_calibration_pearson(mnblist, llist)


#%% MULTIPLE SORTED PLOT RESULTS COMBINED ON ONE PANEL
def sorted_vs_predicted_muitple(fldr_list, clr_list):
    # fixed variables
    num = 100
    axmin, axmax = 0.25, 1.15

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    # axes_lbl_size = 9
    # xtck_size = 9
    # ytck_size = 9
    # leg_size = 8
    # fnt_size = 9
    #
    # plt.rc('font', size=fnt_size)  # controls default text sizes
    # plt.rc('axes', titlesize=axes_lbl_size)  # fontsize of the axes title
    # plt.rc('axes', labelsize=axes_lbl_size)  # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=xtck_size)  # fontsize of the tick labels
    # plt.rc('ytick', labelsize=ytck_size)  # fontsize of the tick labels
    # plt.rc('legend', fontsize=leg_size)  # legend fontsize
    ax1 = plt.subplot(111)
    format_panels(ax1)

    # load data
    for (fldr, col) in zip(fldr_list, clr_list):
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, num)
        div, pi, pred = np.loadtxt(sort_file).T
        rst = cst_from_fldr(fldr)
        # normalize by pi0
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        pi /= pi0
        pred /= pi0
        # plot predicted vs. observed
        lbl = r'B$\geq${:.1f}'.format(0.001*float(fldr[-3:]))
        plt.plot(pred, pi, marker='o', alpha=0.7, ms=2, lw=0.5, label=lbl,
                 color=col)

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='darkslategray', ls='--', alpha=0.65)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(axmin+0.02, 1.02, 'without linked selection', ha='left',
             va='center')

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.3, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    xtick = np.arange(0.3, 1.2, 0.1)
    plt.xticks(xtick, y=0.02)
    plt.ylim(axmin, axmax)
    plt.xlim(axmin, axmax)
    # plt.legend(loc='lower right', ncol=1, frameon=1,
    #            framealpha=0.75, facecolor='white')
    plt.legend(loc='lower right', frameon=True, framealpha=0.75,
               facecolor='white')

    f_save = final_dir + '/sfigs/sorted_compare_B_thresholds.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


bth_tkn = '000 200 500 600'.split()
flist = ['cadd93_bth{}'.format(b) for b in bth_tkn]
clist = ['darkorange', 'dodgerblue', 'mediumvioletred', 'darkturquoise',
         'fuchsia']
sorted_vs_predicted_muitple(flist, clist)


#%% MULTIPLE OBS *MINUS* PRED PLOT RESULTS COMBINED ON ONE PANEL
def obs_minus_pred(fldr_list, clr_list, lbl_list, sname):
    # fixed variables
    num = 100
    axmin, axmax = 0.45, 1.02
    ymin, ymax = -0.07, 0.22
    # create new plot
    plt.figure(figsize=(6.5, 2.5))
    plt.subplots_adjust(top=0.995, right=0.995, left=0.125, bottom=0.155)
    ax1 = plt.subplot(111)
    format_panels(ax1)

    # load data
    for (fldr, col, lbl) in zip(fldr_list, clr_list, lbl_list):
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, num)
        div, pi, pred = np.loadtxt(sort_file).T
        rst = cst_from_fldr(fldr)
        # normalize by pi0
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        pi /= pi0
        pred /= pi0
        # plot predicted - observed
        # lbl = r'B$\geq${:.1f}'.format(0.001*float(fldr[-3:]))
        plt.plot(pred, pi-pred, marker='o', alpha=0.7, ms=3, lw=0.75,
                 label=lbl, color=col)

    plt.axhline(0, color='darkslategray', ls='--', alpha=0.65)
    plt.text(0.5, -0.01, r'$y_i-B_i=0$', va='top', ha='center',
             color='darkslategray', alpha=0.65)
    plt.ylabel('observed minus predicted\n'+r'($y_i-B_i$)', labelpad=3)
    # ytick = np.arange(0.3, 1.2, 0.1)
    plt.yticks(x=0.01)
    plt.xlabel(r'predicted diversity ($B_i$)', labelpad=3)
    xtick = np.arange(0.3, 1.01, 0.1)
    plt.xticks(xtick, y=0.02)
    # plt.ylim(-0.11, 0.33)
    plt.ylim(ymin, ymax)

    plt.xlim(axmin, axmax)
    # plt.legend(loc='lower right', ncol=1, frameon=1,
    #            framealpha=0.75, facecolor='white')
    plt.legend(loc='upper center', frameon=True, framealpha=0.75,
               facecolor='white', ncol=4)

    f_save = final_dir + '/sfigs/{}.obs_minus_pred.png'.format(sname)
    plt.savefig(f_save, dpi=512)
    plt.close()


#
mnb_all = '119 890 678 590 512 442 378 318 263 165'.split()
bth_all = '200 300 400 450 500 550 600 650 700 800'.split()
mnb_tkn_all = [r'B$\geq${:.2f}'.format(float(b)*0.001) for b in bth_all]
llh_tkn_all = [r'LLH$\geq${:.2f}'.format(float(b)*0.001) for b in bth_all]

# bth_tkn = '000 200 500 600'.split()
# mnb_tkn = '119 512 378'.split()
clist = ['darkorange', 'dodgerblue', 'mediumvioletred', 'darkturquoise',
         'fuchsia', 'orangered', 'cyan', 'purple', 'rosybrown']

#
idx = [0, 4, 6]
# anno = 'ape_cons94'
anno = 'cadd94'
# anno = 'fish_cons94'

bthlist = ['{}_gmask_bth_000'.format(anno)]
bthlist += ['{}_gmask_bth_{}'.format(anno, bth_all[i]) for i in idx]
mnblist = ['{}_gmask_bth_000'.format(anno)]
mnblist += ['{}_gmask_mnb_{}'.format(anno, mnb_all[i]) for i in idx]
llist1 = ['none'] + [llh_tkn_all[i] for i in idx]
llist2 = ['none'] + [mnb_tkn_all[i] for i in idx]
sname1 = anno + '_bth'
sname2 = anno + '_mnb'

# obs_minus_pred(bthlist, clist, llist1, sname1)
obs_minus_pred(mnblist, clist, llist2, sname2)

