__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, cst_from_fldr, chromosome_length, \
    human_autosomes, ChromStruct
from classes.phylotree import parse_exptotsub
from data_processing.data_tools import cg_hypermutable, cpg_mask
from figures.common_functions import format_panels, get_bbins, \
    get_telomere_dict, get_centromere_dict, predict_loess
from matplotlib.patches import Rectangle


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


def get_loess_cons(fldr, span):
    fcon = final_dir + '/{}/sorted.cons.n2000.txt'.format(fldr)
    pred, n1, cn1, n2, cn2 = np.loadtxt(fcon)[:2000].T
    wts = np.ones(shape=len(pred))
    cnlo1 = predict_loess(pred, cn1, wts, span, pred)
    cnlo2 = predict_loess(pred, cn2, wts, span, pred)

    return cnlo1, cnlo2


# LOESS FUNCTION FOR SORTMORE DATA
def get_loess_line(fldr, span, return_points=False, load_con=True):
    # load results in 2000 bins for LOESS plots
    fdir = final_dir + '/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_il_n2000.txt'
    div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    # load sorted conserved separately (done using 0.05cM radius of sites)
    if load_con:
        fcon = final_dir + '/{}/sorted.cons.n2000.txt'.format(fldr)
        cn = np.loadtxt(fcon)[:2000,2]

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    points = [pi, gc, cn, cm, il]
    wts = np.ones(shape=len(pred))
    loess_lines = []
    for a in points:
        lo_line = predict_loess(pred, a, wts, span, pred)
        loess_lines.append(lo_line)

    if return_points:
        return pred, loess_lines, points
    else:
        return pred, loess_lines


# COUNT NEUTRAL SITES IN CG HYPERMUTABLE REGIONS PER
def get_cg_masks():
    """get masks of neutral sites in CG hypermutable for each chrom"""
    cg_dict = cg_hypermutable()
    nstr = root_dir + '/data/nmsk/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'
    cg_masks = {}
    cst = ChromStruct('chr1')
    for ch in human_autosomes:
        cst.chrom = ch
        nmsk = np.load(nstr.format(ch))['neutmask']
        c_arr = np.zeros(shape=cst.chlen)
        for (si, sj) in cg_dict[ch]:
            c_arr[si:sj] = 1
        # turn non-neutral sites to 0s
        c_arr[~nmsk] = 0
        cg_masks[ch] = c_arr

    return cg_masks


def get_cg_counts(bbins, cg_masks):
    """get CG hypermutable neutral site counts in 2000 B bins"""
    cg_counts = []
    for (i, b) in enumerate(bbins):
        n_sites = 0
        for (c, start, end, num) in b:
            ch = 'chr{}'.format(int(c))
            start, end = int(start), int(end)
            n_cg = cg_masks[ch][start:end].sum()
            if n_cg > num:
                print i, n_cg, ch, start, end, num
                break
            n_sites += n_cg
        cg_counts.append(n_sites)

    return cg_counts


# GET METHYLATION LEVELS IN B BINS
def get_meth_arrays():
    """get arrays of CpG methylation levels at neutral sites for each chrom"""
    # make empty chrom-length arrays
    meth_dict = {}
    for ch in human_autosomes:
        c_arr = np.zeros(shape=chromosome_length(ch))
        meth_dict[ch] = c_arr

    # fill in methylation levels from Ipsita's bed file across chroms
    # "Methylation levels (measured as percentage methylated)" - Gao et al. 2019
    f_cpgbed = root_dir + '/data/coords/cpgs_for_David.bed'
    with open(f_cpgbed, 'r') as f:
        for line in f:
            ch, start, end, score = line.split()
            start, end = int(start), int(end)
            score = float(score)
            # add scores to chromosome score arrays for every covered site
            if ch in meth_dict:
                meth_dict[ch][start:end] = score

    # flip all methylation scores to 0 at non-neutral sites
    nstr = root_dir + '/data/nmsk/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'
    cst = ChromStruct('chr1')
    for ch in human_autosomes:
        cst.chrom = ch
        # get neutral mask
        nmsk = np.load(nstr.format(ch))['neutmask']
        # turn non-neutral sites to 0s
        meth_dict[ch][~nmsk] = 0

    return meth_dict


def get_meth_sums(bbins, meth_dict):
    """get sum of methylation levels in each of 2000 B bins"""
    methylation_sums = []
    for (i, b) in enumerate(bbins):
        meth_sum = 0
        cpgs_sum = 0
        for (c, start, end, num) in b:
            ch = 'chr{}'.format(int(c))
            start, end = int(start), int(end)
            meth_sum += meth_dict[ch][start:end].sum()
        #     # get the window of meth_dict scores
        #     mscores = meth_dict[ch][start:end]
        #     # get the indices where there are non-zero scores
        #     midx = (mscores > 0)
        #
        #     # sum the CG sites with scores
        #     if np.any(midx):
        #         msum = np.sum(mscores[midx])
        #         # msum = meth_dict[ch][start:end].sum()
        #         # if n_cg > num:
        #         #     print i, n_cg, ch, start, end, num
        #         #     break
        #         meth_sum += msum
        #         cpgs_sum += np.sum(midx)
        # methylation_sums.append(meth_sum / cpgs_sum)
        methylation_sums.append(meth_sum)
    return methylation_sums


# LOAD SUBSTITUTION RATES BY TYPE
def get_sub_rates(fldr, nbins=100, divide_all=False):
    # set paths to data
    rst = cst_from_fldr(fldr)
    bbin_dir = root_dir + '/data/phast/bbins/{}bins'.format(nbins)
    res = []
    for bbin in range(nbins):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        f_name = '{}/{}.bbin{}.exptotsub'.format(bbin_dir, rst.tkn, bbin)

        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m
        # adjust total for number of branches
        for i in range(4):
            mat[i,i] /= 14.0

        # get totals for AT and GC content
        at_tot = (mat[0, :].sum() + mat[3, :].sum())
        cg_tot = (mat[1, :].sum() + mat[2, :].sum())

        # divide by AT/CG by ALL base types
        if divide_all:
            all_tot = at_tot + cg_tot
            at_tot = cg_tot = all_tot

        # get A>C/T>G
        actg = (mat[0, 1] + mat[3, 2]) / at_tot
        # get A>G/T>C
        agtc = (mat[0, 2] + mat[3, 1]) / at_tot
        # get A>T/T>A
        atta = (mat[0, 3] + mat[3, 0]) / at_tot

        # get C>A/G>T
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


# GET SNP POLYMORPHISM BY TYPE
def get_snp_poly(fldr, nbins):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_il_n{}.txt'.format(nbins)
    f_ssnp = fdir + 'sort_snptype_fix_anc_n{}.txt'.format(nbins)

    snp = np.loadtxt(f_ssnp)
    div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

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


#%% ORIGINAL 4 PANEL PLOT FOR TEMPLATES
def original_four_panel(fldr, span=0.1):
    """create loess smoothed plots for conservation and recombination rates"""

    # get 100 points data
    fdir = final_dir + '/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
    div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # get loess line
    prlo, lolist = get_loess_line(fldr, span)
    pilo, gclo, cnlo, cmlo, illo = lolist

    # create new plot
    plt.figure(figsize=(6.5, 6.5))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.1, bottom=0.08,
                        hspace=0.15, wspace=0.3)

    # plot standard predicted/observed on top
    ax1 = plt.subplot(221)
    format_panels(ax1)
    axmin, axmax = 0.55, 1.2
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.97, 0.48
    # plot predicted vs. observed
    plt.plot(pred, pi, marker='o', ms=5, color='darkorange', lw=0,
             alpha=0.5)
    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='darkslategray', ls='--', alpha=0.65)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    # plot LOESS line
    plt.plot(prlo, pilo, lw=2, color='orangered')
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    # plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    xtick = np.arange(0.5, 1.2, 0.1)
    plt.xticks(xtick, y=0.02, color='none')
    plt.ylim(axmin, axmax)
    plt.xlim(axmin, 1.02)
    # plt.legend(loc='lower center', bbox_to_anchor=(0.3, 0.05))
    plt.text(0.75, 0.81, r'$y=x$', rotation=45, ha='center', va='center',
             color='darkslategray', alpha=0.65)
    plt.text(tx1, ty1, 'A', transform=plt.gcf().transFigure)

    # plot recombination rates
    ax2 = plt.subplot(222)
    format_panels(ax2)
    plt.plot(pred, np.log10(cm * 1e6), marker='o', ms=5, lw=0, alpha=0.5,
             color='deepskyblue')
    plt.plot(prlo, np.log10(cmlo*1e6), color='dodgerblue', lw=2)

    plt.xticks(xtick, color='none')
    plt.ylabel('recombination rate (log10 cM/Mb)', labelpad=3)
    # plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylim(-1.2, 1.5)
    plt.xlim(axmin, 1.02)
    plt.text(tx2, ty1, 'B', transform=plt.gcf().transFigure)

    # plot conservation levels
    ax3 = plt.subplot(223)
    format_panels(ax3)
    plt.plot(pred, cn/1000.0, marker='o', ms=5, lw=0, alpha=0.5,
             color='purple')
    plt.plot(prlo, cnlo/1000.0, color='indigo')

    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel('conservation score', labelpad=3)
    plt.xlim(axmin, 1.02)
    # plt.legend(loc='upper left')
    plt.text(tx1, ty2, 'C', transform=plt.gcf().transFigure)

    # plot GC fraction
    ax4 = plt.subplot(224)
    format_panels(ax4)
    plt.plot(pred, gc, marker='o', ms=5, lw=0, alpha=0.5,
             color='forestgreen')
    plt.plot(prlo, gclo, color='darkgreen')

    plt.xticks(xtick)
    plt.ylabel('GC fraction', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.xlim(axmin, 1.02)
    plt.text(tx2, ty2, 'D', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/{}.covars.with.LOESS.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_mnb_378'
original_four_panel(fldr)


#%% OBSERVED VS. PREDICTED INCLUDING CLOSEUP
def observed_vs_predicted(fldr, span=0.1):
    """create loess smoothed plots for conservation and recombination rates"""
    # get 100 points data
    fdir = final_dir + '/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
    div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # get loess line and original points
    prlo, lolist = get_loess_line(fldr, span, load_con=False)
    pilo, gclo, cnlo, cmlo, illo = lolist

    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.07, bottom=0.12,
                        hspace=0.15, wspace=0.3)

    # plot standard predicted/observed on top
    ax1 = plt.subplot(121)
    format_panels(ax1)
    axmin, axmax = 0.55, 1.2
    tx1, tx2 = 0.08, 0.6
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
    # plt.xlabel(r'predicted $\pi/\pi_0$ (B)', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.52, axmax)
    plt.xlim(axmin, 1.02)
    # solve y=x rotation
    adlen = axmax - 0.52
    oplen = adlen * (1.02 - axmin) / (axmax - 0.52)
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')
    plt.text(tx1, ty1, 'a', transform=plt.gcf().transFigure, fontweight='bold')

    # add red box
    inset_lw = 0.75
    pct = 0.9
    rect = Rectangle(xy=(pct, pct), width=0.1, height=0.285,
                     linewidth=inset_lw, edgecolor='crimson',
                     facecolor='none', zorder=5)
    ax1.add_patch(rect)

    # plot closeup
    ax2 = plt.subplot(122)
    format_panels(ax2)
    # get loess line and original points
    prlo, lolist, polist = get_loess_line(fldr, span/2, return_points=True,
                                          load_con = False)
    pilo, gclo, cnlo, cmlo, illo = lolist
    pipo = polist[0]  # pi value points

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='k', ls='--', alpha=1)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    # plot scatter points
    plt.plot(prlo[1400:], pipo[1400:], marker='o', ms=5, color='darkorange',
             lw=0, alpha=0.5)
    # plot loess line
    plt.plot(prlo[1400:], pilo[1400:], lw=2, color='orangered')

    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.801, 1.249)
    plt.xlim(0.903, 0.997)
    plt.text(tx2, ty1, 'b', transform=plt.gcf().transFigure, fontweight='bold')

    f_save = final_dir + '/sfigs/{}.observed.vs.pred.LOESS.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat'
# fldr = 'cadd94_gmask_mnb_378'
observed_vs_predicted(fldr)


#%% RECOMBINATION RATE AND CONSERVATION LEVELS SORTED
def recombination_and_conservation(fldr, span=0.1):
    """create loess smoothed plots for conservation and recombination rates"""

    # get 100 points data
    fdir = final_dir + '/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
    div, pi, pred, gc, _, cm, il = np.loadtxt(f_sort).T

    # load sorted conserved separately (done using 0.05cM radius of sites)
    fcon = final_dir + '/{}/sorted.cons.n100.txt'.format(fldr)
    cn1, cn2 = np.loadtxt(fcon)[:,(2, 4)].T
    cnlo1, cnlo2 = get_loess_cons(fldr, span)

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # get loess line
    prlo, lolist = get_loess_line(fldr, span, load_con=False)
    pilo, gclo, cnlo, cmlo, illo = lolist

    # create new plot
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(top=0.98, right=0.99, left=0.102, bottom=0.16,
                        hspace=0.15, wspace=0.35)

    # plot standard predicted/observed on top
    axmin, axmax = 0.55, 1.2
    tx1, tx2, tx3 = 0.104, 0.427, 0.752
    ty1, ty2 = 0.92, 0.48
    xtick = np.arange(0.5, 1.2, 0.1)
    xtcky = 0.04
    ytckx = 0.04
    xlpad = 1
    ylpad = 1

    # 1. plot recombination rates
    ax1 = plt.subplot(131)
    format_panels(ax1)
    plt.plot(pred, np.log10(cm * 1e6), marker='o', ms=5, lw=0, alpha=0.5,
             color='deepskyblue')
    plt.plot(prlo, np.log10(cmlo * 1e6), color='dodgerblue', lw=2)

    plt.xticks(xtick)
    # plt.ylabel('recombination rate (cM/Mb in log scale)', labelpad=3)
    plt.ylabel('recombination rate\n' + r'$\mathrm{(log_{10}(cM/Mb))}$',
               labelpad=ylpad)

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=xlpad)
    plt.ylim(-1.2, 1.4)
    # ytckvalues = [-1, -0.5, 0, 0.5, 1]
    # ytckstrings = [r'$10^{%.1f}$' %y for y in ytckvalues]
    # plt.yticks(ytckvalues, ytckstrings, x=0.03)
    plt.yticks(x=ytckx)
    plt.xticks(y=xtcky)
    plt.xlim(axmin, 1.02)
    plt.text(tx1, ty1, 'a', transform=plt.gcf().transFigure, fontweight='bold')

    # 2. plot conservation levels (in 50kb windows)
    ax2 = plt.subplot(132)
    format_panels(ax2)

    plt.plot(pred, cn2/1000.0, marker='o', ms=5, lw=0, alpha=0.5,
             color='darkturquoise')
    plt.plot(prlo, cnlo2/1000.0, color='cadetblue')
    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick, y=xtcky)
    plt.yticks(x=ytckx)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=xlpad)
    plt.ylabel('functional density per bp', labelpad=ylpad)
    # plt.ylabel('mean phastCons score\n(50 kb radius per site)',
    #            labelpad=ylpad)
    plt.xlim(axmin, 1.02)
    # plt.legend(loc='upper left')
    plt.text(tx2, ty1, 'b', transform=plt.gcf().transFigure, fontweight='bold')

    # 3. plot conservation levels (in 0.05cM windows)
    ax3 = plt.subplot(133)
    format_panels(ax3)
    plt.plot(pred, cn1 / 1000.0, marker='o', ms=5, lw=0, alpha=0.5,
             color='purple')
    plt.plot(prlo, cnlo1 / 1000.0, color='indigo')

    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=xlpad)
    plt.ylabel('functional density per cM', labelpad=ylpad)
    # plt.ylabel('mean phastCons score\n(0.05 cM radius per site)',
    #            labelpad=ylpad)

    # plt.ylabel('conservation score (in 0.05 cM radius)', labelpad=3)
    plt.xlim(axmin, 1.02)
    # plt.yticks(x=ytckx)
    plt.xticks(y=xtcky)
    ytcks = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14]
    plt.yticks(ytcks, x=ytckx)
    plt.ylim(0.025, 0.145)
    # plt.legend(loc='upper left')
    plt.text(tx3, ty1, 'c', transform=plt.gcf().transFigure, fontweight='bold')

    f_save = final_dir + '/sfigs/{}.recomb.and.cons.LOESS.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()

fldr = 'cadd94_gmask_v1.6_without_bstat'
# fldr = 'cadd94_gmask_mnb_378'
recombination_and_conservation(fldr)


#%% GC CONTENT AND GC IN CPG ISLANDS
def gc_content(fldr, span=0.1, show_diff=False):
    """create loess smoothed plots for conservation and recombination rates"""

    # get 100 points data
    fdir = final_dir + '/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
    div, pi, pred, gc, _, cm, il = np.loadtxt(f_sort).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # get loess line
    prlo, lolist = get_loess_line(fldr, span, load_con=False)
    pilo, gclo, cnlo, cmlo, illo = lolist

    # tick and letter positions
    xtick = np.arange(0.5, 1.2, 0.1)
    axmin, axmax = 0.55, 1.2
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.94, 0.48
    xlpad = 1
    ylpad = 1
    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.1, bottom=0.13,
                        hspace=0.15, wspace=0.3)

    # plot GC fraction
    ax1 = plt.subplot(121)
    format_panels(ax1)
    plt.plot(pred, gc*100.0, marker='o', ms=5, lw=0, alpha=0.5,
             color='forestgreen')
    plt.plot(prlo, gclo*100.0, color='darkgreen', label='all GC')

    # show difference when content in islands is removed
    if show_diff:
        plt.plot(pred, gc-il, marker='o', ms=5, lw=0, alpha=0.5,
                 color='lime')
        plt.plot(prlo, gclo-illo, color='limegreen',
                 label='remove GC in islands')
        plt.legend()

    plt.xticks(xtick)
    plt.ylabel('GC content (%)', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xlim(axmin, 1.02)
    plt.text(tx1, ty1, 'a', transform=plt.gcf().transFigure, fontweight='bold')

    # plot GC fraction in CpG islands
    ax2 = plt.subplot(122)
    format_panels(ax2)
    # NOTE: multiple by 2 for counting both the C and the G, since counts
    # are for the Cs only. multiply by 100 to get in units %
    plt.plot(pred, 200*il, marker='o', ms=5, lw=0, alpha=0.5,
             color='mediumorchid')
    plt.plot(prlo, 200*illo, color='mediumpurple')
    print np.mean(200*il)

    plt.xticks(xtick)
    plt.ylabel('proportion of CpGs in CpG islands (%)', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xlim(axmin, 1.02)

    plt.text(tx2, ty1, 'b', transform=plt.gcf().transFigure, fontweight='bold')
    # f_label = '{}.GC.with.LOESS'.format(fldr)
    # f_label = '{}.ancestral.GC.with.LOESS'.format(fldr)
    f_label = '{}.ancestral.GC.with.LOESS.CpGs.in.islands'.format(fldr)

    if show_diff:
        f_label += '.showdiff'
    f_save = final_dir + '/sfigs/{}.png'.format(f_label)
    plt.savefig(f_save, dpi=512)
    plt.close()


# fldr = 'cadd94_gmask_mnb_378'
fldr = 'cadd94_gmask_v1.6_without_bstat'
gc_content(fldr, show_diff=False)


#%% SPATIAL DISTRIBUTION OF NEUTRAL SITES SORTED BY B
def get_spatial_distributions(bins):
    """
    Spatial distributions as described in figure caption:

    In (A), we measure the distance to telomeres on the same chromosome arm
    In (B), we show the distribution of relative chromosomal position, i.e.,
    the distance to the centromere relative to the distance between
    centromere and telomeres on that chromosomal arm, for all putatively
    neutral sites and for those in the top 1% of B-values; relative distances
    on the shorter arm are to the left
    """
    # chrom lookup dictionaries for telomere/centromere coordinates
    centro_dict = get_centromere_dict()  # centromere dictionary
    telo_dict = get_telomere_dict()  # telomere dictionary

    # lists for data across B bins
    mean_dist = []  # mean distances and stddev
    rel_dist_sh = []  # relative distances on short arm
    rel_dist_ln = []  # relative distances on long arm
    n_site_sh = []  # counts on short arm
    n_site_ln = []  # counts on long arm

    # empty lists for B~1 data
    b1_dsh = []  # short arm dists
    b1_dln = []  # long arm dists
    b1_nsh = []  # short arm counts
    b1_nln = []  # long arm coutns

    # counter for sites that were inside centromere
    inside_ctm = 0

    # process data from each B bin
    for (n, b) in enumerate(bins):
        # set data lists for current B bin
        d_abs = []
        d_srel = []
        d_lrel = []
        n_site = []
        n_sh = []
        n_ln = []
        for (c, start, end, num) in b:
            # set current chrom name
            ch = 'chr{}'.format(int(c))
            # collect the number of neutral sites in the current segment
            n_site.append(num)
            # get centromere start/end for the current chrom
            ctm_start, ctm_end = centro_dict[ch]
            # get telomere points for current chrom
            tlm_start, tlm_end = telo_dict[ch]
            # get the length of the short arm of the chrom
            if (ctm_start - tlm_start) > (tlm_end - ctm_end):
                short_arm = tlm_end - ctm_end
            else:
                short_arm = ctm_start - tlm_start
            # get mean point for current segment
            mean_point = start + ((end - start) / 2.0)

            # if point on 5' arm, get distance to 5' telomere and arm length
            if mean_point < ctm_start:
                dist = mean_point - tlm_start
                arm_len = ctm_start - tlm_start
            # if point on 3' arm, get distance to 3' telomere and arm length
            elif mean_point > ctm_end:
                dist = tlm_end - mean_point
                arm_len = tlm_end - ctm_end
            # if inside centromere, place on closest side and record #
            else:
                inside_ctm += 1
                ctm_mid = 1.0 * (ctm_end - ctm_start)
                if mean_point > ctm_mid:
                    dist = tlm_end - mean_point
                    arm_len = tlm_end - ctm_mid
                else:
                    dist = mean_point - tlm_start
                    arm_len = ctm_mid - tlm_start

            # record the absolute distance to the telomere
            d_abs.append(dist)

            # calculate the relative distance from the centromere
            rel_dist = float(arm_len - dist) / arm_len
            assert 0 <= rel_dist <= 1

            # record relative distance, sorted by short/long chrom arms
            if arm_len == short_arm:
                d_srel.append(rel_dist)
                n_sh.append(num)
            else:
                d_lrel.append(rel_dist)
                n_ln.append(num)

        # gather all data for the B bin
        rel_dist_sh.extend(d_srel)
        rel_dist_ln.extend(d_lrel)
        n_site_sh.extend(n_sh)
        n_site_ln.extend(n_ln)

        # get the mean and (weighted) stddev of the absolute distances
        m_dist = np.average(d_abs, weights=n_site)
        v_dist = np.average((np.array(d_abs)-m_dist)**2, weights=n_site)
        std_dist = np.sqrt(v_dist)
        mean_dist.append((m_dist, std_dist))

        # get the data for the final B~1 bin
        if n == 99:
            b1_dsh = d_srel
            b1_dln = d_lrel
            b1_nsh = n_sh
            b1_nln = n_ln

    # print the number of aberrant mid-centromere sites there were
    print 'there were {} points found within centromeres'.format(inside_ctm)

    # get histogram of sites for each spatial bin
    dbins = np.arange(0, 1.01, 0.01)
    all_rel_sh = np.histogram(rel_dist_sh, bins=dbins, weights=n_site_sh)[0]
    all_rel_ln = np.histogram(rel_dist_ln, bins=dbins, weights=n_site_ln)[0]
    b1_sh = np.histogram(b1_dsh, bins=dbins, weights=b1_nsh)[0]
    b1_ln = np.histogram(b1_dln, bins=dbins, weights=b1_nln)[0]

    # return histogram counts and mean distances + stddev
    return dbins, all_rel_sh, all_rel_ln, b1_sh, b1_ln, mean_dist


# tkn = 'cadd94_gmask'
tkn = 'cadd94_gmask_v1.6_without_bstat'
flbls = ['dbins', 'allrelshort', 'allrellong', 'b1short', 'b1long', 'meandist']
pfiles = [final_dir + '/{}/{}.txt'.format(tkn, fl) for fl in flbls]
if not all(os.path.isfile(pf) for pf in pfiles):
    print 'making relative position files for {}.'.format(tkn)
    bins = get_bbins(tkn)
    db, rsh, rln, b1sh, b1ln, mdist = get_spatial_distributions(bins)
    # save these values so they don't need to be rerun
    data = [db, rsh, rln, b1sh, b1ln, mdist]
    for i in range(6):
        np.savetxt(pfiles[i], data[i])
else:
    print 'loading relative position files for {}.'.format(tkn)
    db, rsh, rln, b1sh, b1ln, mdist = [np.loadtxt(pf) for pf in pfiles]


#%% PLOT RELATIVE POSITIONS OF NEUTRAL SITES ON CHROM ARMS
def plot_relative_chrom_position(all_rel_sh, all_rel_ln, b1_sh, b1_ln):
    """plot relative chrom position of B stratified sites, with B~1 separate"""
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.995, right=0.995, bottom=0.11, left=0.16)
    ax = plt.subplot(111)
    format_panels(ax)
    xleft = np.arange(-1, 0, 0.01)
    xright = np.arange(0.01, 1.01, 0.01)
    plt.step(xleft, 1.0 * all_rel_sh[::-1] / np.sum(all_rel_sh),
             color='dodgerblue', alpha=0.75,
             label='all putatively neutral sites')
    plt.step(xright, 1.0 * all_rel_ln / np.sum(all_rel_ln), color='dodgerblue',
             alpha=0.75)
    plt.step(xleft, 1.0 * b1_sh[::-1] / np.sum(b1_sh), color='forestgreen',
             alpha=0.75)
    plt.step(xright, 1.0 * b1_ln / np.sum(b1_ln), color='forestgreen',
             alpha=0.75, label=r'100th percentile of predicted $B$')

    plt.xticks(y=0.02)
    plt.xlabel('relative chromosomal position', labelpad=2)
    plt.yticks(x=0.02)
    plt.ylabel('proportion of neutral sites', labelpad=2)
    plt.ylim(-0.002, 0.142)
    plt.text(0.18, 0.94, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.legend(loc='upper center')

    f_save = final_dir + '/sfigs/relative.chromosomal.positions.png'
    plt.savefig(f_save, dpi=516)
    plt.close()


plot_relative_chrom_position(rsh, rln, b1sh, b1ln)


#%% PLOT THE DISTANCE TO THE TELOMERES FOR NEUTRAL SITES STRATIFIED BY B BIN
def plot_distance_to_telomere(fldr, mean_dist):
    f_sort = final_dir + '/{}/basic_sort_n100.txt'.format(fldr)
    div, pi, pred = np.loadtxt(f_sort).T
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pred /= pi0

    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.995, right=0.995, bottom=0.115, left=0.145)
    mn, st = np.array(mean_dist).T
    ax1 = plt.subplot(111)
    format_panels(ax1)
    plt.errorbar(pred, mn, yerr=st/2, ecolor='deepskyblue', fmt='o',
                 label=r'mean $\pm$ std', color='dodgerblue')
    plt.xticks(y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xlim(0.55, 1.02)
    plt.yticks(x=0.02)
    plt.ylim(0, 7.9e7)
    # plt.ylabel(r'distance to telomere ($\times10^7$bp)', labelpad=2)
    plt.ylabel('average distance to telomere\n(in units of 10 Mb)', labelpad=2)

    plt.legend(loc='upper center')
    plt.text(0.16, 0.94, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/{}.distance.to.telomeres.by.B.png'.format(fldr)
    plt.savefig(f_save, dpi=516)
    plt.close()


plot_distance_to_telomere('cadd94_gmask_v1.6_without_bstat', mdist)


#%% PLOT METHYLATION % AND % IN C>G HYPERMUTABLE
def plot_methdist_and_cghyperdist(fldr):
    # get prediction data
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pred_list = []
    for nbins in [100, 2000]:
        f_sort = final_dir + '/{}/'.format(
            fldr) + 'sort_gc_cm_cn_il_n{}.txt'.format(nbins)
        pred = np.loadtxt(f_sort)[:, 2]
        pred /= pi0
        pred_list.append(pred)
    pr100, pr2000 = pred_list

    # load or generate and save CG data
    cgcount_file = final_dir + '/{}/cg_count_2000.txt'.format(fldr)
    if os.path.isfile(cgcount_file):
        print 'loading CG data.'
        cgsums = np.loadtxt(cgcount_file)
    else:
        print 'generating and saving CG data.'
        cg_msk = get_cg_masks()
        bbin = get_bbins(fldr, 2000)
        cgsums = get_cg_counts(bbin, cg_msk)
        np.savetxt(cgcount_file, cgsums)
    # get CG hypermutable LOESS and 100 bins points
    cglo = predict_loess(pr2000, cgsums, None, 0.1, pr2000)
    cg100 = []
    for i in range(0, 2000, 20):
        cg_avg = np.sum(cgsums[i:i + 20])
        cg100.append(cg_avg)
    cg100 = np.array(cg100)

    # load or generate and save methylation data
    methsum_file = final_dir + '/{}/meth_count_2000.txt'.format(fldr)
    if os.path.isfile(methsum_file):
        print 'loading methylation data.'
        msums = np.loadtxt(methsum_file)
    else:
        print 'generating and saving methylation data.'
        mdict = get_meth_arrays()
        bbin = get_bbins(fldr, 2000)
        msums = get_meth_sums(bbin, mdict)
        np.savetxt(methsum_file, msums)
    # get methlyation LOESS and 100 bins points
    mthlo = predict_loess(pr2000, msums, None, 0.1, pr2000)
    mth100 = []
    for i in range(0, 2000, 20):
        m = np.sum(msums[i:i + 20])
        mth100.append(m)
    mth100 = np.array(mth100)

    # tick and letter positions
    xtick = np.arange(0.5, 1.2, 0.1)
    axmin, axmax = 0.55, 1.2
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.94, 0.48
    xlpad = 1
    ylpad = 1

    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.1, bottom=0.13,
                        hspace=0.15, wspace=0.3)

    neut_per2000 = 653e6 / 2000.0
    neut_per100 = 653e6 / 100.0
    ax1 = plt.subplot(121)
    format_panels(ax1)
    plt.plot(pr100, 100*mth100/neut_per100, color='brown', alpha=0.5, lw=0,
             marker='o', ms=5)
    plt.plot(pr2000, 100*mthlo/neut_per2000, color='maroon', lw=2)

    plt.xticks(y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.yticks(x=0.02)
    plt.xlim(axmin, 1.02)

    plt.ylabel('methylated CpG content (%)', labelpad=2)
    # plt.ylim(-0.01, 0.24)
    plt.text(tx1, ty1, 'c', transform=plt.gcf().transFigure, fontweight='bold')

    # plt.legend(loc='upper center')

    ax2 = plt.subplot(122)
    format_panels(ax2)
    plt.plot(pr100, 100*cg100/neut_per100, color='LightSeaGreen', alpha=0.5, lw=0,
             marker='o', ms=5)
    plt.plot(pr2000, 100*cglo/neut_per2000, color='CadetBlue', lw=2)

    plt.xticks(y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xlim(axmin, 1.02)

    plt.yticks(x=0.02)
    plt.ylabel('neutral sites in C>G hypermutable (%)', labelpad=2)
    # plt.ylim(-0.01, 0.24)
    plt.text(tx2, ty1, 'd', transform=plt.gcf().transFigure, fontweight='bold')


    f_save = final_dir + '/sfigs/mthCpG.and.CG.hypermut.png'
    plt.savefig(f_save, dpi=516)
    plt.close()


# plot_methylation_distribution(pr100, msums100, pr2000, methlo)
# plot_cg_hypermutable_distribution(pr100, cg_cnt100, pr2000, cglo)
fldr = 'cadd94_gmask_v1.6_without_bstat'
plot_methdist_and_cghyperdist(fldr)


#%% GET 100 AND 2000 BIN DATA FOR STANDARD DATA
fldr = 'cadd94_gmask_v1.6_without_bstat'
# fldr = 'cadd94_gmask_mnb_378'
n100_poly = get_snp_poly(fldr, 100)
n2000_poly = get_snp_poly(fldr, 2000)
n100_subs = get_sub_rates(fldr)
n2000_subs = get_sub_rates(fldr, 2000)


# SMOOTH 2000 SCALE SUBSTITUTION RATES WITH LOESS
pr2000 = n2000_poly[1]
wts = np.ones(shape=len(pr2000))
span = 0.1
dv_loess = []
for d in n2000_subs[1]:
    lo = predict_loess(pr2000, d, wts, span, pr2000)
    dv_loess.append(lo)


# SMOOTH 2000 SCALE SNPS NORMALIZED BY SUBSTITUTION RATES
pi_loess1 = []
pi_loess2 = []
pi_point = []
for i in range(6):
    p = n2000_poly[0][i]
    d = n2000_subs[1][i]
    lo1 = predict_loess(pr2000, p/d, None, 0.1, pr2000)
    lo2 = predict_loess(pr2000[1400:], p[1400:]/d[1400:], None, 0.2,
                        pr2000[1400:])
    pi_loess1.append(lo1)
    pi_loess2.append(lo2)
    pi_point.append(p[1400:]/d[1400:])


#%% GET 100 AND 2000 BIN DATA FOR CG FILTERED DATA
fldr_cg = 'cadd94_gmask_filter_CG'
# fldr_cg = 'cadd94_gmask_v1.6_without_bstat_filter_CG'
n100_poly_cg = get_snp_poly(fldr_cg, 100)
n2000_poly_cg = get_snp_poly(fldr_cg, 2000)
n100_subs_cg = get_sub_rates(fldr_cg)
n2000_subs_cg = get_sub_rates(fldr_cg, 2000)

# # make out of the smaller bins
# n100_subs_cg = []
# for i in range(0, 2000, 20):
#     n100_subs_cg.append([n2000_subs_cg[1][j][i:i+20].mean() for j in range(6)])
# n100_subs_cg = np.array(n100_subs_cg).T

# SMOOTH 2000 SCALE SNPS NORMALIZED BY SUBSTITUTION RATES
pr2000_cg = n2000_poly_cg[1]
pi_loess1_cg = []
pi_loess2_cg = []
pi_point_cg = []
for i in range(6):
    p = n2000_poly_cg[0][i]
    d = n2000_subs_cg[1][i]
    lo1 = predict_loess(pr2000_cg, p/d, None, 0.1, pr2000_cg)
    lo2 = predict_loess(pr2000_cg[1400:], p[1400:]/d[1400:], None, 0.2,
                        pr2000_cg[1400:])
    pi_loess1_cg.append(lo1)
    pi_loess2_cg.append(lo2)
    pi_point_cg.append(p[1400:]/d[1400:])


#%% RESHUFFLE BINNED DATA BY % C>G HYPERMUTABLE
fldr = 'cadd94_gmask_v1.6_without_bstat'
neut_per100 = 653e6 / 100.0
# get prediction data
rst = cst_from_fldr(fldr)
pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
pred_list = []
for nbins in [100, 2000]:
    f_sort = final_dir + '/{}/'.format(
        fldr) + 'sort_gc_cm_cn_il_n{}.txt'.format(nbins)
    pred = np.loadtxt(f_sort)[:, 2]
    pred /= pi0
    pred_list.append(pred)
pr100, pr2000 = pred_list

cgcount_file = final_dir + '/{}/cg_count_2000.txt'.format(fldr)
if os.path.isfile(cgcount_file):
    print 'loading CG data.'
    cgsums = np.loadtxt(cgcount_file)
else:
    print 'generating and saving CG data.'
    cg_msk = get_cg_masks()
    bbin = get_bbins(fldr, 2000)
    cgsums = get_cg_counts(bbin, cg_msk)
    np.savetxt(cgcount_file, cgsums)
# get CG hypermutable LOESS and 100 bins points
cglo = predict_loess(pr2000, cgsums, None, 0.1, pr2000)
cg100 = []
for i in range(0, 2000, 20):
    cg_avg = np.sum(cgsums[i:i + 20])
    cg100.append(cg_avg)
cg100 = np.array(cg100)
cg100 /= neut_per100
# calculate the fraction of 20 bins to keep (in 2000 bin data)
frac = (20 * (1-cg100)).astype(int)
# get neut poly and subs data for 2000 bins
n2000_poly = get_snp_poly(fldr, 2000)  # ancpi, pred, snp
n2000_subs = get_sub_rates(fldr, 2000)  # res, dv

#
api, prd = n2000_poly[:2]
dv = n2000_subs[1]
shpi, shpr, shdv = [], [], []
idx = 0
for i in range(0, 2000, 20):
    ii = i+frac[idx]
    shpr.append(prd[i:ii])
    curpi, curdv = [], []
    for j in range(6):
        curpi.append(api[j][i:ii])
        curdv.append(dv[j][i:ii])
    shpi.append(curpi)
    shdv.append(curdv)
    idx += 1
#
newpi = []
newdv = []
newpr = np.concatenate(shpr)
for i in range(6):
    temppi = []
    tempdv = []
    for j in range(100):
        temppi.extend(shpi[j][i])
        tempdv.extend(shdv[j][i])
    newpi.append(temppi)
    newdv.append(tempdv)

stepsize = int(len(newpr) / 100)
shpi, shdv, shpr = [], [], []
for i in range(0, len(newpr), stepsize):
    shpr.append(newpr[i:i+stepsize].mean())
    curpi, curdv = [], []
    for j in range(6):
        curpi.append(np.mean(newpi[j][i:i + stepsize]))
        curdv.append(np.mean(newdv[j][i:i + stepsize]))
    shpi.append(curpi)
    shdv.append(curdv)


shpi = np.array(shpi).T
shdv = np.array(shdv).T


#%% PLOT COMPARING SHUFFLED AND REGULAR CG FILTERED
plt.figure(figsize=(9.75, 6.5))
plt.subplots_adjust(top=0.96, right=0.99, hspace=0.2, bottom=0.06,
                    left=0.07)
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
         'rosybrown']
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
# plt.title('Ancestral AT originating mutations')
x_limits = (0.58, 1.01)
ylimits = (0.47, 1.41)
pred = n100_poly_cg[1]
# AT SNPS
for i in [0, 1, 2]:
    ax1 = plt.subplot(2, 3, i+1)
    format_panels(ax1)
    # plot regular C>G filtered
    pi = n100_poly_cg[0][i] / n100_subs_cg[1][i]
    plt.plot(pred, pi, marker='s', ms=5, lw=0, alpha=1, markerfacecolor='None',
             markeredgecolor=clist[i], markeredgewidth=0.9,
             label='C>G filtered')
    # plot reshuffled
    pi = shpi[i] / shdv[i]
    plt.plot(shpr, pi, marker='o', ms=5, lw=0, alpha=0.5, color=clist[i],
             label='random filtered')
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.ylim(0.55, 1.15)
    plt.xlim(*x_limits)
    plt.ylim(*ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white', title=labs[i])

# CG SNPS
for i in [3, 4, 5]:
    ax1 = plt.subplot(2, 3, i+1)
    format_panels(ax1)
    # plot regular C>G filtered
    pi = n100_poly_cg[0][i] / n100_subs_cg[1][i]
    plt.plot(pred, pi, marker='s', ms=5, lw=0, alpha=1, markerfacecolor='None',
             markeredgecolor=clist[i], markeredgewidth=0.9,
             label='C>G filtered')
    # plot reshuffled
    pi = shpi[i] / shdv[i]
    plt.plot(shpr, pi, marker='o', ms=5, lw=0, alpha=0.5, color=clist[i],
             label='random filtered')
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.ylim(0.55, 1.15)
    plt.xlim(*x_limits)
    plt.ylim(*ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white', title=labs[i])

f_save = final_dir + '/sfigs/CG-hyper-vs-shuffled_AT_CG_SNPs.png'
plt.savefig(f_save, dpi=512)
plt.close()


#%% PLOT COMPARING SHUFFLED AND REGULAR CG FILTERED (2)
plt.figure(figsize=(6, 6))
plt.subplots_adjust(top=0.96, right=0.99, hspace=0.2, bottom=0.12,
                    left=0.13)
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
         'rosybrown']

x_limits = (0.59, 1.01)
ylimits = (0.47, 1.15)

# plot standard predicted/observed on top
ax1 = plt.subplot(111)
format_panels(ax1)
axmin, axmax = 0.55, 1.2
tx1, tx2 = 0.08, 0.6
ty1, ty2 = 0.94, 0.48
xtick = np.arange(0.5, 1.2, 0.1)

# get 100 points data
fldr = 'cadd94_gmask_v1.6_without_bstat'
fdir = final_dir + '/{}/'.format(fldr)
f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T
rst = cst_from_fldr(fldr)
pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
pi /= pi0
pred /= pi0

# plot y=x line
plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
         color='k', ls='--')
# plot horizontal line at y=1
plt.axhline(y=1, color='k', alpha=0.8, ls='-')
# plot no filtering
plt.plot(pred, pi, marker='^', ms=5, color='dodgerblue', lw=0,
         alpha=0.5, label='no filtering')
# plot regular C>G filtered
pred = n100_poly_cg[1]
pi = np.average(n100_poly_cg[0], axis=0) / np.average(n100_subs_cg[1], axis=0)
plt.plot(pred, pi, marker='s', ms=5, lw=0, alpha=1, markerfacecolor='None',
         markeredgecolor='darkorange', markeredgewidth=0.9,
         label='C>G filtered')
# plot reshuffled
pi = np.average(shpi, axis=0) / np.average(shdv, axis=0)
plt.plot(shpr, pi, marker='o', ms=5, lw=0, alpha=0.5, color='purple',
         label='random filtered')
plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
# plt.xticks()
# plt.yticks()
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
# plt.ylim(0.55, 1.15)
plt.xlim(*x_limits)
plt.ylim(*ylimits)
plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white')

f_save = final_dir + '/sfigs/CG-hyper-vs-shuffled-combine-SNPs.png'
plt.savefig(f_save, dpi=512)
plt.close()


#%% AT AND CG ORIGINATING SUBSTITUTIONS
def plot_at_cg_substitutions(fldr, pred, dv, predlo, dvlo):
    # set labels and colors for each class of substitution
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown']
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.07, bottom=0.12,
                        hspace=0.15, wspace=0.3)

    # set positions for A-B labels and axes limits
    tx1, tx2 = 0.08, 0.6
    ty1, ty2 = 0.94, 0.48
    x_limits = (0.58, 1.02)
    sub_ylimits = (0.83, 1.45)

    # AT SUBSTITUTIONS
    ax3 = plt.subplot(121)
    format_panels(ax3)
    # plot 100 bins points first
    for i in [0, 1, 2]:
        plt.plot(pred, dv[i], marker='o', ms=5, lw=0, color='white')
        plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.4, color=clist[i])

    # plot 2000 bins loess smoothed line next
    for i in [0, 1, 2]:
        plt.plot(predlo, dvlo[i], label=labs[i], lw=2, color=clist[i])

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'relative substitution rate ($D/\bar{D}$)', labelpad=2)
    plt.xlim(*x_limits)
    plt.ylim(*sub_ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white')
    plt.text(tx1, ty1, 'a', transform=plt.gcf().transFigure, fontweight='bold')

    # CG SUBSTITUTIONS
    ax4 = plt.subplot(122)
    format_panels(ax4)
    # plot 100 bins points first
    for i in [3, 4, 5]:
        plt.plot(pred, dv[i], marker='o', ms=5, lw=0, color='white')
        plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.4, color=clist[i])

    # plot 2000 bins loess smoothed line next
    for i in [3, 4, 5]:
        plt.plot(predlo, dvlo[i], label=labs[i], lw=2, color=clist[i])

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xlim(0.55, 1.02)
    plt.ylabel(r'relative substitution rate ($D/\bar{D}$)', labelpad=2)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.xlim(*x_limits)
    plt.ylim(*sub_ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white')
    plt.text(tx2, ty1, 'b', transform=plt.gcf().transFigure, fontweight='bold')

    f_save = final_dir + '/sfigs/{}_AT_CG_substitutions.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


plot_at_cg_substitutions(fldr, n100_poly[1], n100_subs[1], pr2000, dv_loess)


#%% AT AND CG ORIGINATING POLYMORPHISM
def plot_at_cg_snps(fldr, pred, ancpi, dv, predlo, pilo1, pilo2, pipt):
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown']
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    plt.figure(figsize=(6.5, 6.5))
    plt.subplots_adjust(top=0.96, right=0.99, hspace=0.2, bottom=0.06,
                        left=0.07)

    tx1, tx2 = 0.075, 0.577
    ty1, ty2 = 0.935, 0.44
    x_limits = (0.58, 1.01)
    if 'CG' in fldr:
        ylimits = (0.47, 1.41)
        ylim_closeup = (0.6, 1.55)

    else:
        ylimits = (0.47, 1.65)
        ylim_closeup = (0.6, 1.65)

    lo_alpha = 1
    # AT SNPS
    ax1 = plt.subplot(221)
    format_panels(ax1)
    plt.title('Ancestral AT originating mutations')

    # y=x line
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    # plot 100 bins points first
    for i in [0, 1, 2]:
        pi = ancpi[i] / dv[i]
        plt.plot(pred, pi, marker='o', ms=5, lw=0, color='white')
        plt.plot(pred, pi, marker='o', ms=5, lw=0, alpha=0.4, color=clist[i])

    # plot LOESS line next
    for i in [0, 1, 2]:
        plt.plot(predlo, pilo1[i], label=labs[i], lw=2, color=clist[i],
                 alpha=lo_alpha)

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.ylim(0.55, 1.15)
    plt.xlim(*x_limits)
    plt.ylim(*ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white')
    plt.text(tx1, ty1, 'a', transform=plt.gcf().transFigure, fontweight='bold')
    # add red box
    inset_lw = 0.75
    pct = 0.9
    rect = Rectangle(xy=(pct, pct), width=0.1, height=0.2,
                     linewidth=inset_lw, edgecolor='crimson',
                     facecolor='none', zorder=10)
    ax1.add_patch(rect)

    # CG SNPS
    ax2 = plt.subplot(222)
    format_panels(ax2)
    plt.title('Ancestral CG originating mutations')

    # y=x line
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    # plot 100 bins points first
    for i in [3, 4, 5]:
        pi = ancpi[i] / dv[i]
        plt.plot(pred, pi, marker='o', ms=5, lw=0, color='white')
        plt.plot(pred, pi, marker='o', ms=5, lw=0, alpha=0.4, color=clist[i])

    # plot LOESS line next
    for i in [3, 4, 5]:
        plt.plot(predlo, pilo1[i], label=labs[i], lw=2, color=clist[i],
                 alpha=lo_alpha)

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.ylim(0.55, 1.15)
    plt.xlim(*x_limits)
    plt.ylim(*ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75,
               facecolor='white')
    plt.text(tx2, ty1, 'b', transform=plt.gcf().transFigure, fontweight='bold')
    # add red box
    inset_lw = 0.75
    pct = 0.9
    if 'CG' in fldr:
        hgt = 0.3
    else:
        hgt = 0.7
    rect = Rectangle(xy=(pct, pct), width=0.1, height=hgt,
                     linewidth=inset_lw, edgecolor='crimson',
                     facecolor='none', zorder=10)
    ax2.add_patch(rect)

    # AT CLOSEUP
    ax3 = plt.subplot(223)
    format_panels(ax3)
    # y=x line
    plt.plot([0.9, 1], [0.9, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    # plot points
    for i in [0, 1, 2]:
        plt.plot(predlo[1400:], pipt[i], marker='o', ms=3, lw=0,
                 color='white')
        plt.plot(predlo[1400:], pipt[i], marker='o', ms=3, lw=0,
                 alpha=0.4, color=clist[i])
    # plot LOESS line next
    for i in [0, 1, 2]:
        plt.plot(predlo[1400:], pilo2[i], label=labs[i], lw=2,
                 color=clist[i], alpha=lo_alpha)

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.xlim(0.903, 0.997)
    plt.ylim(*ylim_closeup)
    # plt.legend(loc='upper center', frameon=1, framealpha=0.75, facecolor='white')
    plt.text(tx1, ty2, 'c', transform=plt.gcf().transFigure, fontweight='bold')

    # CG CLOSEUP
    ax4 = plt.subplot(224)
    format_panels(ax4)
    # y=x line
    plt.plot([0.9, 1], [0.9, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    #plot points
    for i in [3, 4, 5]:
        plt.plot(predlo[1400:], pipt[i], marker='o', ms=3, lw=0,
                 color='white')
        plt.plot(predlo[1400:], pipt[i], marker='o', ms=3, lw=0,
                 alpha=0.4, color=clist[i])
    # plot LOESS line next
    for i in [3, 4, 5]:
        plt.plot(predlo[1400:], pilo2[i], label=labs[i], lw=2,
                 color=clist[i], alpha=lo_alpha)

    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.xlim(0.903, 0.997)
    plt.ylim(*ylim_closeup)
    # plt.legend(loc='upper center', frameon=1, framealpha=0.75, facecolor='white')
    plt.text(tx2, ty2, 'd', transform=plt.gcf().transFigure, fontweight='bold')

    f_save = final_dir + '/sfigs/{}_AT_CG_SNPs.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


# plot_at_cg_snps(fldr, n100_poly[1], n100_poly[0], n100_subs[1], pr2000,
#                 pi_loess1, pi_loess2, pi_point)
plot_at_cg_snps(fldr_cg, n100_poly_cg[1], n100_poly_cg[0], n100_subs_cg[1],
                pr2000_cg, pi_loess1_cg, pi_loess2_cg, pi_point_cg)


#%% UN-SCALED HETEROZYGOSITY AND SUBSTITUTIONS BY MUTATION TYPE
fldr = 'cadd94_gmask_v1.6_without_bstat'
span = 0.1
sub_loess = []
n2000_total_subs = get_sub_rates(fldr, 2000, divide_all=True)
n100_total_subs = get_sub_rates(fldr, 100, divide_all=True)
for s in n2000_total_subs[0].T:
    lo = predict_loess(pr2000, s, None, span, pr2000)
    sub_loess.append(lo)

snp_loess = []
n100_total_poly = get_snp_poly(fldr, 100)
n2000_total_poly = get_snp_poly(fldr, 2000)
for p in n2000_total_poly[2].T:
    lo = predict_loess(pr2000, p, None, span, pr2000)
    snp_loess.append(lo)


# UN-SCALED HETEROZYGOSITY AND SUBSTITUTIONS BY MUTATION TYPE
def unscaled_pi_and_div(fldr, pred, snp, sub, prlo, snplo, sublo):
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown']
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.92, right=0.99, hspace=0.1, bottom=0.11,
                        left=0.07, wspace=0.25)
    txt_y = 0.865
    plt.suptitle('Breakdown of SNPs and substitutions by type')
    # PLOT SNPS
    ax1 = plt.subplot(121)
    format_panels(ax1)
    # plt.title('Heterozygosity by SNP')
    # plot the points for diversity for each SNP type
    for (i,s) in enumerate(snp[:,:6].T):
        plt.plot(pred, s*1e3, marker='o', ms=3, lw=0, alpha=0.5,
                 color=clist[i])
    # plot the LOESS lines from 2000 bins for each SNP type
    for i in range(6):
        plt.plot(prlo, snplo[i]*1e3, label=labs[i], color=clist[i], lw=1)

    # xtick = np.arange(0.6, 1.01, 0.1)
    # plt.xticks(xtick)
    plt.xticks(y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.yticks(x=0.03)
    plt.ylabel(r'contribution to diversity level $\times 10^3$', labelpad=2)
    plt.xlim(0.58, 1.02)
    plt.ylim(0.0, 0.87)
    plt.legend(loc='upper center', ncol=2, handletextpad=0.2, borderpad=0.2,
               borderaxespad=0.3, frameon=1, handlelength=1,
               framealpha=0.75, facecolor='white', columnspacing=0.4,
               labelspacing=0.4)
    plt.text(0.08, txt_y, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')

    # PLOT SUBS
    ax2 = plt.subplot(122)
    format_panels(ax2)
    # plt.title('Substitution rates')

    # plot the points for divergence for each substitution type
    for (i, r) in enumerate(sub[:, :6].T):
        plt.plot(pred, r, marker='o', ms=3, lw=0, alpha=0.5, color=clist[i])

    # plot the LOESS lines from 2000 bins for each sub type
    for i in range(6):
        plt.plot(prlo, sublo[i], label=labs[i], color=clist[i], lw=1)

    # xtick = np.arange(0.5, 1.01, 0.1)
    # plt.xticks(xtick)
    plt.xticks(y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.yticks(x=0.03)
    plt.ylabel('substitutions per site', labelpad=2)
    # plt.ylabel('contribution to divergence', labelpad=2)
    plt.xlim(0.58, 1.02)
    plt.ylim(0, 0.072)
    # plt.legend(loc='upper left', ncol=2)
    plt.text(0.592, txt_y, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/{}.unscaled.SNPs.and.subs.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


unscaled_pi_and_div(fldr, n100_total_poly[1], n100_total_poly[2],
                    n100_total_subs[0], pr2000, snp_loess, sub_loess)


#%% PLOT R^2 AND OBS/PRED FOR CADD VS. CADD FILTER CG
def plot_rsq_and_sortpred(flist, clist, llist):
    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.07, bottom=0.12,
                        hspace=0.15, wspace=0.3)

    # plot standard predicted/observed on top
    ax1 = plt.subplot(121)
    format_panels(ax1)
    axmin, axmax = 0.55, 1.2
    tx1, tx2 = 0.08, 0.6
    ty1, ty2 = 0.94, 0.48
    xtick = np.arange(0.5, 1.2, 0.1)

    pi_list, pr_list, pilo_list, prlo_list = [], [], [], []
    for fl in flist:
        # get 100 points data
        fdir = final_dir + '/{}/'.format(fl)
        f_sort = fdir + 'sort_gc_cm_cn_il_n100.txt'
        div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

        rst = cst_from_fldr(fl)
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        pi /= pi0
        pred /= pi0
        pi_list.append(pi)
        pr_list.append(pred)

        # get loess line and original points
        span = 0.1
        prlo, lolist = get_loess_line(fl, span, load_con=False)
        pilo = lolist[0]
        pilo_list.append(pilo)
        prlo_list.append(prlo)

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='k', ls='--', alpha=1)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    # plot predicted vs. observed
    for i in range(len(flist)):
        plt.plot(pr_list[i], pi_list[i], marker='o', ms=5, color=clist[i][0],
                 lw=0, alpha=0.5)

    # plot LOESS line
    for i in range(len(flist)):
        plt.plot(prlo_list[i], pilo_list[i], lw=2, color=clist[i][1],
                 label=llist[i])

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
    plt.text(tx1, ty1, 'b', transform=plt.gcf().transFigure, fontweight='bold')
    # plt.legend(loc='lower right')

    # PLOT R^2 RESULT
    ax2 = plt.subplot(122)
    format_panels(ax2)
    for i in range(len(flist)):
        f_rsq = final_dir + '/{}/rsq.log'.format(flist[i])
        w, r = np.loadtxt(f_rsq)[:16].T
        plt.plot(np.log10(w), r, marker='o', ms=5, color=clist[i][0], lw=0,
                 alpha=0.75)
    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x % 1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.02)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.text(tx2, ty1, 'c', transform=plt.gcf().transFigure, fontweight='bold')

    f_save = final_dir + '/sfigs/CADD_filter_CG_sort_rsq.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


# flist = ['cadd94_gmask_mnb_378', 'cadd94_gmask_filter_CG']
flist = ['cadd94_gmask_v1.6_without_bstat',
         'cadd94_gmask_v1.6_without_bstat_filter_CG']
clist = [['darkorange','orangered'], ['purple', 'indigo']]
llist = ['CADD', 'CADD without\nC>G hypermutable']
plot_rsq_and_sortpred(flist, clist, llist)


#%% COMPARE PARAMS FROM THE TWO RUNS
def compare_params(flist, llist, clist):
    # create lists for data to plot
    dfe_list = []
    pmf_list = []
    pi0_list = []

    for fldr in flist:
        # initalize ChromStruct from final composite from folder
        cst = cst_from_fldr(fldr)

        # scale upmf by u0
        upmf = [u/u0 for u in cst.uvec]
        pmf_list.append(upmf)
        # get pi0 based on tau
        mpi0 = cst.params[-1] / cst.fixed.tau_init
        pi0_list.append(mpi0)
        # use representative DFE for all values
        dfe_list.append(cst.bdfe[0])

    # plot formatting parameters
    xi = np.arange(len(dfe_list[0]))
    xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]
    n = sum(len(pm) for pm in pmf_list)
    w = 0.8 / n
    s = -w * n / 2.0
    ncol = 1

    # create the plot
    plt.figure(figsize=(6.5, 2))
    plt.subplots_adjust(left=0.07, wspace=1.2, right=0.995, bottom=0.17, top=0.91)
    alpha_val = 0.4
    shaded_col = 'white'
    plt.rc('hatch', color='k', linewidth=0.25)
    hatch_pattern = '\\'*6
    # DFE plot
    ax1 = plt.subplot(1, 6, (1, 4))
    format_panels(ax1)

    plt.title('a', fontweight='bold', loc='left', y=0.975)
    plt.title('(i)', fontweight='bold', loc='center', y=0.975)
    ymax = 0

    # solid bar results
    istop = len(pmf_list) / 2
    leg1lines = []
    for (i, df) in enumerate(pmf_list):
        assert len(df) == 1
        lbl = llist[i]
        col = clist[i]
        df = df[0]
        ymax = max(df.max(), ymax)
        plt.bar(xi + s, df, w, label=lbl, color=col, align='edge')
        s += w

    plt.legend(loc='upper left', ncol=ncol, borderaxespad=0.3,borderpad=0.2,
                      handlelength=1.2, labelspacing=0.25)

    # # hatched bar results
    # leg2lines = []
    # for (i, df) in enumerate(pmf_list[istop:]):
    #     assert len(df) == 1
    #     lbl = llist[i]
    #     col = clist[i]
    #     df = df[0]
    #     ymax = max(df.max(), ymax)
    #     li = plt.bar(xi + s, df, w, color=col, label=lbl, align='edge',
    #             hatch=hatch_pattern)[0]
    #     leg2lines.append(li)
    #     s += w
    # plt.legend(leg2lines, llist[istop:], loc='upper right', ncol=ncol,
    #            title='CADD 6%', borderaxespad=0.3,borderpad=0.2,
    #            handlelength=1.2, labelspacing=0.25)
    # plt.gca().add_artist(leg1)

    # make small ticks between selection coefs
    for xx in xi[:-1]:
        plt.axvline(xx + 0.5, 0, 0.2, lw=1, ls='--', color='k')

    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.02)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    plt.xlabel('deleterious fitness effect', labelpad=3)

    # plot udel
    ax2 = plt.subplot(165)
    format_panels(ax2)

    plt.title('(ii)', loc='center', y=0.975, fontweight='bold')
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    s = -w / 2.0
    li = 0
    for (i, df) in enumerate(pmf_list):
        for d in df:
            udel = sum(d)
            plt.bar(0 + s, udel, w, color=clist[li])
            li += 1
            s += w
    plt.xticks([])
    plt.yticks(x=0.15)

    # plot pi0
    ax3 = plt.subplot(166)
    format_panels(ax3)

    plt.title('(iii)', loc='center', y=0.975, fontweight='bold')
    # plt.ylabel('average reduction\nin heterozygosity', labelpad=3)
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=3)

    plt.yticks(x=0.15)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        # if i%2:
        plt.bar(0 + s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    f_save = final_dir + '/sfigs/updateJune2021.CG.filter.params.png'
    plt.savefig(f_save, dpi=512)


flist = ['cadd94_gmask_v1.6_without_bstat',
         'cadd94_gmask_v1.6_without_bstat_filter_CG']
clist = ['darkorange', 'purple']
llist = ['CADD', 'CADD without\nC>G hypermutable']
compare_params(flist, llist, clist)


#%% LOAD ARCHAIC INTROGRESSION LEVELS SORTED BY B
def get_archaic(fldr, pop, num):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_arch = fdir + 'predsort_archaic_{}_n{}.txt'.format(pop, num)
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pred, ai = np.loadtxt(f_arch).T
    pred /= pi0

    return pred, ai


def archaic_introgression_plot(fldr, pop, ttl, letter):
    pr100, ai100 = get_archaic(fldr, pop, 100)
    pr2000, ai2000 = get_archaic(fldr, pop, 2000)
    ailo = predict_loess(pr2000, ai2000, None, span, pr2000)
    fig_dir = root_dir + '/result/final_files/sfigs/'
    f_save = fig_dir + '{}_archaic_{}.png'.format(fldr, pop)

    plt.figure(figsize=(3.25,3.25))
    plt.subplots_adjust(top=0.93, right=0.99, hspace=0.1, bottom=0.14,
                        left=0.25)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title(ttl, y=0.98)
    plt.plot(pr100, ai100, marker='o', ms=5, lw=0, alpha=0.5,
             color='forestgreen')

    # plot LOESS curve
    plt.plot(pr2000, ailo, color='darkgreen', label='all GC')


    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylabel('estimated proportion of NE ancestry', labelpad=3)
    plt.xlim(0.55, 1.02)
    plt.text(0.26, 0.88, letter, transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat'
pops = ['CEU', 'CHBS']
titles = ['European (CEU)', 'East-Asian (CHB/CHS)']
letters = ['a', 'b']
for i in [0,1]:
    archaic_introgression_plot(fldr, pops[i], titles[i], letters[i])
#%%
