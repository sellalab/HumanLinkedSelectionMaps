__author__ = 'davidmurphy'

import os
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir
from figures.other_code.summary_slide import cst_from_fldr
from figures.common_functions import format_panels
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


#%% CHROM 1 FIGURES
def chr1_diversity(show_pred=False):
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.07, bottom=0.155, right=0.995, top=0.995)
    fldr = 'cadd94_gmask_mnb_378'
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    chlist = [f for f in os.listdir(fdir) if ('chr1' in f) and
             (f.endswith('.txt')) and 'TMRCA' not in f]
    f_name = fdir + chlist[0]
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan

    ax = plt.subplot(111)
    format_panels(ax)
    plt.plot(xi, obs, label='observed', color='darkslategray', lw=2)
    if show_pred:
        plt.plot(xi, prd, color='darkorange', lw=1.5, alpha=0.8,
                 label='predicted')
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    plt.ylim(0, 1.9)
    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x / 25) % 2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])

    # figure letter
    if show_pred:
        tkn = 'with_pred'
        loc = 'lower center'
        plt.legend(loc=loc, ncol=2, frameon=1,
                   framealpha=0.75, facecolor='white', handlelength=0.8,
                   borderaxespad=0.3, columnspacing=0.8)
    else:
        tkn = 'without_pred'
    f_save = final_dir + '/sfigs/chr1.{}.png'.format(tkn)
    plt.savefig(f_save, dpi=512)
    plt.close()


chr1_diversity()
chr1_diversity(True)


#%%  COLLATED PLOT SLIDES
def collate_plots(show_pred=False, show_mcv=False):
    flist = ['cadd94_gmask_mnb_378', 'mcvicker']

    # plt.figure(figsize=(6.5, 2.16))
    plt.figure(figsize=(4.32, 2.16))
    plt.subplots_adjust(left=0.12, bottom=0.18, right=0.995, top=0.995,
                        wspace=0.3)

    # COLLATED PLOT
    # ax1 = plt.subplot(1,3,(1,2))
    ax1 = plt.subplot(111)
    format_panels(ax1)
    f_our = final_dir + '/{}/nonsyn.collate.5.00e-03width.npy'.format(flist[0])
    bins, div, pi, pr, cnts = np.load(f_our).T
    obs = pi / cnts
    prd = pr / cnts
    y_max = obs.mean()
    obs /= y_max
    prd /= y_max
    bins += 2.5e-3

    f_mcv = root_dir + '/result/collate/YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy'
    mbins, mdiv, mpi, mpr, mcnts = np.load(f_mcv).T
    mobs = mpi / mcnts
    mprd = mpr / mcnts
    my_max = mobs.mean()
    mobs /= my_max
    mprd /= my_max
    mbins += 2.5e-3

    plt.plot(bins, obs, label='observed', color='darkslategray', lw=2)
    if show_pred:
        plt.plot(bins, prd, label='our map', color='darkorange', alpha=0.8,
                 lw=1.5)
    if show_mcv:
        plt.plot(mbins, mprd, label='McVicker map', color='mediumpurple',
                 alpha=0.8,
                 lw=1.5)

    plt.xlabel('distance to nearest NS substitution (cM)', labelpad=3)
    plt.xticks(y=0.03)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi}$)', labelpad=3)
    ytck = np.arange(0.7, 1.2, 0.1)
    plt.yticks(ytck, x=0.02)
    plt.ylim(0.68, 1.15)
    if (show_mcv or show_pred):
        plt.legend(loc='lower right')
    tkn = 'collated'
    if show_pred:
        tkn += '.pred'
    if show_mcv:
        tkn += '.mcv'
    f_save = final_dir + '/sfigs/{}.png'.format(tkn)
    plt.savefig(f_save, dpi=512)
    plt.close()

    # # R^2 PLOT
    # dall = 'fish_cons94_gmask_mnb_378'
    # f_mcv = root_dir + '/result/final_files/mcvicker/rsq.log'
    # f_cad = root_dir + '/result/final_files/{}/rsq.log'.format(dall)
    # w, mcv = np.loadtxt(f_mcv)[:16].T
    # ape = np.loadtxt(f_cad)[:16, 1]
    #
    # w = np.log10(w)
    # # plt.figure(figsize=(3.25, 2.5))
    # ax2 = plt.subplot(133)
    # format_panels(ax2)
    # plt.plot(w, ape, label='our map', marker='o', lw=0,
    #          color='darkorange', ms=5)
    # plt.plot(w, mcv, label='McVicker', marker='o', lw=0, color='mediumpurple',
    #          ms=5)
    # # plt.plot(w, ex, label='our method (exon only)', marker='o', lw=0,
    # #          color='fuchsia')
    #
    # plt.xlabel('window size (bp)', labelpad=1)
    # xtck = [4, 4.5, 5, 5.5, 6]
    # xstr = [r'$10^{%.1f}$' % x if not x%1 else '' for x in xtck]
    # plt.xticks(xtck, xstr, y=0.04)
    # plt.ylabel('variance explained ' + r'$(\mathrm{R^2})$', labelpad=1)
    # plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.04)
    # plt.ylim(0.01, 0.65)
    # # plt.legend(prop=dict(size=9), loc='upper left')
    # plt.text(0.073, 0.92, 'B', transform=plt.gcf().transFigure)
    # plt.text(0.742, 0.92, 'C', transform=plt.gcf().transFigure)

    # f_save = root_dir + '/result/final_files/sfigs/fig_20.rsq.mcvicker.png'
    # plt.savefig(f_save, dpi=256)
    # plt.close()

    # plt.text(0.01, 0.93, fletters[0], transform=plt.gcf().transFigure)
    # plt.text(0.68, 0.93, fletters[1], transform=plt.gcf().transFigure)


collate_plots(0, 0)
collate_plots(0, 1)
collate_plots(1, 0)
collate_plots(1, 1)
#%% R^2 PLOT
def rsq_plot():
    plt.figure(figsize=(2, 2))
    plt.subplots_adjust(left=0.24, bottom=0.2, top=0.99, right=0.97)
    ax = plt.subplot(111)
    format_panels(ax)
    fldr_1 = 'cadd94_gmask_mnb_378'
    fldr_2 = 'fish_cons94_gmask_exonic'
    clist = ['darkorange', 'dodgerblue']
    llist = ['all', 'exonic']
    # load R^2 data for all conserved and exonic conserved (ape)
    rsq_list = []
    flist = [fldr_1, fldr_2]
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(fldr_1))[:16, 0])
    for fl in flist:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    for i in xrange(1):
        rs = rsq_list[i]
        lab = llist[i]
        col = clist[i]
        plt.plot(w, rs, label=lab, color=col, marker='o', lw=0, ms=4, alpha=0.8)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x % 1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.02)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $\mathrm{(R^2)}$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    # plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
    #            borderaxespad=0.3, ncol=1, frameon=1,
    #            framealpha=0.75, facecolor='white', columnspacing=0.1,
    #            labelspacing=0.25)
    # figure letters
    # plt.text(0.075, 0.915, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.7475, 0.915, 'B', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/rsqplot.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


rsq_plot()


#%% R^2 PLOT 2
def rsq_plot_2(flist, llist, clist, sname):
    plt.figure(figsize=(2, 2))
    plt.subplots_adjust(left=0.24, bottom=0.2, top=0.99, right=0.97)
    ax = plt.subplot(111)
    format_panels(ax)
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(flist[0]))[:16, 0])
    for fl in flist:
        fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    for i in xrange(len(flist)):
        rs = rsq_list[i]
        lab = llist[i]
        col = clist[i]
        plt.plot(w, rs, label=lab, color=col, marker='o', lw=0, ms=4, alpha=0.8)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x % 1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.02)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $\mathrm{(R^2)}$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.legend(loc='upper left', handletextpad=-0.5, borderpad=0.1,
               borderaxespad=0.0, ncol=1,
               labelspacing=0.1, prop=dict(size=8))

    f_save = final_dir + '/sfigs/rsq_compare_annos_{}.png'.format(sname)
    plt.savefig(f_save, dpi=512)
    plt.close()


# flist = ['fish_cons94_gmask_mnb_378', 'huvec_chromhmm_gmask',
#          'cadd94_gmask_mnb_378', 'genic_gmask']
# llist = ['conservation', 'ChromHMM', 'CADD', 'genic' ]
# clist = ['darkorange', 'darkturquoise', 'dodgerblue', 'mediumpurple']
# flist = ['fish_cons94_gmask_mnb_378', 'YRI_nonsyn_s1']
# llist = ['background selection', 'selective sweeps']
# clist = ['darkorange', 'darkturquoise']
# flist = ['fish_cons94_gmask_mnb_378', 'mcvicker']
# llist = ['our map', 'McVicker' ]
# clist = ['darkorange', 'mediumpurple']
# flist = ['fish_cons94_gmask_mnb_378', 'fish_cons94_gmask_exonic',
#          'fish_cons94_gmask_nonexonic']
# llist = ['both', 'exonic', 'nonexonic']
# clist = ['darkorange', 'dodgerblue', 'purple']

flist = ['cadd94_gmask_mnb_378', 'genic_plus_ccre', 'genic']
llist = ['CADD', 'coding/cCRE', 'genic']
clist = ['darkorange', 'dodgerblue', 'purple']
rsq_plot_2(flist, llist, clist, 'CADD_codingcCRE_genic')
#%% UDEL SIMPLE
def udel_simple():
    ulow, uhigh = 1.29, 1.51
    mcvicker = 7.4
    cons = 0.988
    cadd = 1

    ui = [mcvicker, cons, cadd]
    xi = np.arange(len(ui))
    cols = ['mediumpurple', 'orange', 'dodgerblue']
    xticks = ['McVicker', 'conservation', 'CADD']

    plt.figure(figsize=(3.25, 2.16))
    plt.subplots_adjust(left=0.1, bottom=0.09, right=0.995, top=0.9)
    ax = plt.subplot(111)
    format_panels(ax)
    ax.grid(lw=0)
    plt.bar(xi, ui, color=cols)
    plt.axhline(ulow, color='darkslategray', ls='-', lw=0.5)
    plt.axhline(uhigh, color='darkslategray', ls='-', lw=0.5)
    fill_lbl = r'$\mu_{tot}$ range'
    plt.fill_between([-1, 3], ulow, uhigh, color='gray', alpha=0.5,
                     label=fill_lbl)

    plt.title('Deleterious mutation rates', y=0.97)
    plt.xlim(-0.5, 2.5)
    plt.xticks(xi, xticks, y=0.04)
    plt.yticks(np.arange(1, 8), x=0.02)
    plt.ylim(0, 7.8)
    plt.ylabel(r'$\mu_{del}\ (\times10^8)$', labelpad=2)
    plt.legend()

    f_save = final_dir+ '/sfigs/udel_simple.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


udel_simple()
#%%