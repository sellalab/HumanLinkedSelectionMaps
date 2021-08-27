__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, cst_from_fldr, chromosome_length
from classes.phylotree import parse_exptotsub
from data_processing.functions import rsquared_function
from figures.common_functions import format_panels, get_bbins, \
    get_telomere_dict, get_centromere_dict

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


#%% LOAD SUBSTITUTION RATES BY TYPE (100 BINS)
def get_sub_rates(fldr, nbins=100):
    rst = cst_from_fldr(fldr)
    bbin_dir = root_dir + '/data/phast/bbins/{}bins'.format(nbins)

    res = []
    for bbin in xrange(nbins):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        f_name = '{}/{}.bbin{}.exptotsub'.format(bbin_dir, rst.tkn, bbin)

        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m
        # adjust total for number of branches
        for i in xrange(4):
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


# GET SUBSTITUTION RATES BY TYPE (TOP X OF 2000 BINS)
def get_sub_rates_2(fldr):
    rst = cst_from_fldr(fldr)
    bbin_dir = root_dir + '/data/phast/bbins'
    res = []
    msk = []
    for bbin in xrange(1500, 2000):
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
        for i in xrange(4):
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


# GET SNP POLYMORPHISM BY TYPE
def get_snp_poly(fldr, nbins=100):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_n{}.txt'.format(nbins)
    f_ssnp = fdir + 'sort_snptype_fix_anc_n{}.txt'.format(nbins)

    snp = np.loadtxt(f_ssnp)
    div, pi, pred = np.loadtxt(f_sort)[:,:3].T
    # if fldr == 'cadd93_extel_rm_CG':
    #     rst = cst_from_fldr('cadd93_bth600')
    # else:
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


#%% SORT PREDICTIONS, CONSERVATION, CMMMB AND CG FRACTION
def sort_data(fldr, nbins=100):
    """sort additional divergence data for B=1 analyses"""
    # load results
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_sort = fdir + 'sort_gc_cm_cn_n{}.txt'.format(nbins)
    div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    # create new plot
    plt.figure(figsize=(6.5, 6.5))
    plt.subplots_adjust(top=0.995, right=0.99, left=0.1, bottom=0.08,
                        hspace=0.15, wspace=0.3)

    # plot standard predicted/observed on top
    ax1 = plt.subplot(221)
    format_panels(ax1)
    axmin, axmax = 0.55, 1.12
    color = 'darkorange'
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.97, 0.48
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
    plt.plot(pred, np.log10(cm * 1e6), marker='o', ms=5, lw=0, alpha=0.75,
             color='deepskyblue', label='CADD 7%')
    plt.xticks(xtick, color='none')
    plt.ylabel('recombination rate (log10 cM/Mb)', labelpad=3)
    # plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylim(-1.2, 1.2)
    plt.xlim(axmin, 1.02)
    plt.text(tx2, ty1, 'B', transform=plt.gcf().transFigure)

    # plot conservation levels
    ax3 = plt.subplot(223)
    format_panels(ax3)
    plt.plot(pred, cn/1000.0, marker='o', ms=5, lw=0, alpha=0.75,
             color='purple')
    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel('conservation score', labelpad=3)
    plt.xlim(axmin, 1.02)
    plt.legend(loc='upper left')
    plt.text(tx1, ty2, 'C', transform=plt.gcf().transFigure)

    # plot GC fraction
    ax4 = plt.subplot(224)
    format_panels(ax4)
    plt.plot(pred, gc, marker='o', ms=5, lw=0, alpha=0.75,
             color='forestgreen')
    plt.xticks(xtick)
    plt.ylabel('GC fraction', labelpad=3)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.xlim(axmin, 1.02)
    plt.text(tx2, ty2, 'D', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/{}.other_sorted_data.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


sort_data('cadd94_gmask_mnb_378', nbins=100)
#%% SPATIAL DISTRIBUTION OF NEUTRAL SITES SORTED BY B
def get_spatial_distributions(bins):
    centro_dict = get_centromere_dict()
    telo_dict = get_telomere_dict()
    b_1_dists = []
    b_1_num = []
    all_b_dists = []
    all_b_num = []
    mean_dist = []
    for (n, b) in enumerate(bins):
        dists = []
        nsite = []
        arm_dists = []
        for (c, start, end, num) in b:
            nsite.append(num)
            ch = 'chr{}'.format(int(c))
            # get centromere start/end for the current chrom
            ctm_start, ctm_end = centro_dict[ch]
            # get telomere points for current chrom
            tlm_start, tlm_end = telo_dict[ch]
            # get mean point for segment
            mean_point = start + ((end - start) / 2.0)
            # check what side of centromere the point falls on
            if mean_point < ctm_start:
                # get distance to centromere
                dist = ctm_start - mean_point
                # calculate the length of the current chrom arm
                arm_len = ctm_start - tlm_start
            elif mean_point > ctm_end:
                # get distance to centromere
                dist = mean_point - ctm_end
                # calculate the length of the current chrom arm
                arm_len = tlm_end - ctm_end
            else:
                print 'ERROR: POINT IS WITHIN CENTROMERE'
                dist, arm_len = 0, 1

            # normalize distance by chromosome arm length
            norm_dist = 1.0 * mean_point / chromosome_length(ch)
            arm_dist = 1.0 * dist / arm_len

            assert norm_dist <= 1
            dists.append(norm_dist)
            arm_dists.append(arm_dist)

        # gather all of the distances into a single list
        all_b_dists.extend(dists)
        all_b_num.extend(nsite)

        # get the data for B=1 sites
        if n == 99:
            b_1_dists = dists
            b_1_num = nsite

        # get the mean and (weighted) stddev of the distances
        mdist = np.average(arm_dists, weights=nsite)
        variance = np.average((np.array(arm_dists)-mdist)**2, weights=nsite)
        stdist = np.sqrt(variance)

        mean_dist.append((mdist, stdist))

    # GET HISTOGRAM OF SITES FOR EACH SPATIAL BIN
    dbins = np.arange(0, 1.01, 0.01)
    all_count = np.histogram(all_b_dists, bins=dbins, weights=all_b_num)[0]
    b1_count = np.histogram(b_1_dists, bins=dbins, weights=b_1_num)[0]

    return dbins, all_count, b1_count, mean_dist


tkn = 'cadd94_gmask'
bins = get_bbins(tkn)
dbins, all_count, b1_count, mean_dist = get_spatial_distributions(bins)


#%% PLOT RELATIVE CHROMOSOMAL POSITION OF ALL SITES VS. B=1 SITES
def plot_b1_outliers(dbins, all_count, b1_count):
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.995, right=0.995, bottom=0.11, left=0.16)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.step(dbins[:-1], 1.0 * all_count / np.sum(all_count),
             color='dodgerblue', alpha=0.75, label='all')
    plt.step(dbins[:-1], 1.0 * b1_count / np.sum(b1_count),
             color='forestgreen', alpha=0.75, label=r'B$\approx$1')
    plt.xticks(y=0.02)
    plt.xlabel('relative chromosomal position', labelpad=2)
    plt.yticks(x=0.02)
    plt.ylabel('proportion of neutral sites', labelpad=2)
    plt.ylim(-0.002, 0.142)
    plt.text(0.18, 0.94, 'B', transform=plt.gcf().transFigure)
    plt.legend(loc='upper center')

    f_save = final_dir + '/sfigs/fig_S28.dist_all_b1.png'
    plt.savefig(f_save, dpi=516)
    plt.close()


plot_b1_outliers(dbins, all_count, b1_count)


#%% PLOT BMAP SORTED RELATIVE POSITIONS OF NEUTRAL SITES W/ STDDEV
def plot_spatial_distribution(fldr, mean_dist):
    f_sort = final_dir + '/{}/basic_sort_n100.txt'.format(fldr)
    div, pi, pred = np.loadtxt(f_sort).T
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pred /= pi0

    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.995, right=0.995, bottom=0.11, left=0.135)
    mn, st = np.array(mean_dist).T
    ax1 = plt.subplot(111)
    format_panels(ax1)
    plt.errorbar(pred, mn, yerr=st/2, ecolor='dodgerblue', fmt='o',
                 label=r'mean distance $\pm$ std')
    plt.xticks(y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    plt.xlim(0.55, 1.02)
    plt.yticks(x=0.02)
    plt.ylabel('normalized distance to centromere', labelpad=2)
    plt.legend(loc='upper center')
    plt.text(0.16, 0.94, 'A', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/{}.spatial.bdist.png'.format(fldr)
    plt.savefig(f_save, dpi=516)
    plt.close()


plot_spatial_distribution('cadd94_gmask_mnb_378', mean_dist)


#%%
# fldr = 'cadd94_gmask_mnb_378'
fldr = 'cadd94_gmask_filter_CG'

res, dv = get_sub_rates(fldr)
ancpi, pred, snp = get_snp_poly(fldr, 100)

#%% AT AND CG ORIGINATING SNPS AND SUBSTITUTIONS
def plot_at_cg_mutation(fldr, pred, ancpi, dv):
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown']
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    plt.figure(figsize=(6.5, 6.5))
    plt.subplots_adjust(top=0.96, right=0.99, hspace=0.1, bottom=0.06,
                        left=0.07)

    tx1, tx2 = 0.075, 0.577
    ty1, ty2 = 0.935, 0.46
    x_limits = (0.58, 1.02)
    sub_ylimits = (0.83, 1.37)
    snp_ylimits = (0.49, 1.22)
    # AT SNPS
    ax1 = plt.subplot(221)
    format_panels(ax1)
    plt.title('AT originating mutations')
    for i in [0, 1, 2]:
        plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.xticks(color='none')
    plt.yticks(x=0.02)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.ylim(0.55, 1.15)
    plt.xlim(*x_limits)
    plt.ylim(*snp_ylimits)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75, facecolor='white')
    plt.text(tx1, ty1, 'A', transform=plt.gcf().transFigure)

    # CG SNPS
    ax2 = plt.subplot(222)
    format_panels(ax2)
    plt.title('CG originating mutations')
    for i in [3, 4, 5]:
        plt.plot(pred, ancpi[i] / dv[i], marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)

    plt.xlim(*x_limits)
    plt.ylim(*snp_ylimits)
    # plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.legend(loc='upper center', frameon=1, framealpha=0.75, facecolor='white')
    plt.text(tx2, ty1, 'B', transform=plt.gcf().transFigure)

    # AT SUBSTITUTIONS
    ax3 = plt.subplot(223)
    format_panels(ax3)
    for i in [0, 1, 2]:
        plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.ylabel(r'$D/\bar{D}$', labelpad=2)
    plt.xlim(*x_limits)
    plt.ylim(*sub_ylimits)
    # plt.legend(loc='upper center', frameon=1, framealpha=0.75, facecolor='white')
    plt.text(tx1, ty2, 'C', transform=plt.gcf().transFigure)

    # CG SUBSTITUTIONS
    ax4 = plt.subplot(224)
    format_panels(ax4)
    for i in [3, 4, 5]:
        plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    plt.xlim(0.55, 1.02)
    # plt.ylabel(r'$D/\bar{D}$', labelpad=2)
    plt.xticks(y=0.02)
    plt.yticks(x=0.02)
    plt.xlim(*x_limits)
    plt.ylim(*sub_ylimits)
    # plt.legend(loc='upper center', frameon=1, framealpha=0.75, facecolor='white')
    plt.text(tx2, ty2, 'D', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/{}_AT_CG_mutations.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


plot_at_cg_mutation(fldr, pred, ancpi, dv)


#%% UN-SCALED HETEROZYGOSITY AND SUBSTITUTIONS BY MUTATION TYPE
def plot_unscaled_at_cg_mutation(fldr, snp, res, pred):
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown']
    labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
    plt.figure(figsize=(6.5, 3.2))
    plt.subplots_adjust(top=0.92, right=0.99, hspace=0.1, bottom=0.11,
                        left=0.07, wspace=0.25)
    txt_y = 0.87

    # PLOT SNPS
    ax1 = plt.subplot(121)
    format_panels(ax1)
    plt.title('Heterozygosity by SNP')
    for (i,s) in enumerate(snp[:,:6].T):
        plt.plot(pred, s*1e3, marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    # xtick = np.arange(0.6, 1.01, 0.1)
    # plt.xticks(xtick)
    plt.xticks(y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    plt.yticks(x=0.03)
    plt.ylabel(r'observed $\pi$ ($\times 10^3$)', labelpad=2)
    plt.xlim(0.58, 1.02)
    plt.ylim(0.02, 1.75)
    plt.legend(loc='upper center', ncol=2,handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1,
               labelspacing=0.25)

    plt.text(0.075, txt_y, 'A', transform=plt.gcf().transFigure)

    # PLOT SUBS
    ax2 = plt.subplot(122)
    format_panels(ax2)
    plt.title('Substitution rates')
    for (i, r) in enumerate(res[:, :6].T):
        plt.plot(pred, r, marker='o', ms=5, lw=0, alpha=0.75,
                 label=labs[i], color=clist[i])
    # xtick = np.arange(0.5, 1.01, 0.1)
    # plt.xticks(xtick)
    plt.xticks(y=0.03)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    plt.yticks(x=0.03)
    plt.ylabel('substitution rate', labelpad=2)
    plt.xlim(0.58, 1.02)
    # plt.legend(loc='upper left', ncol=2)
    plt.text(0.585, txt_y, 'B', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/{}.unscaled.mutations.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


plot_unscaled_at_cg_mutation(fldr, snp, res, pred)


#%% PARAMS FOR CG FILTERED VERSION OF CADD 6%
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
    plt.rc('hatch', color='darkgray', linewidth=0.25)

    # DFE plot
    ax1 = plt.subplot(1, 6, (1, 4))
    format_panels(ax1)

    plt.title('A (i)', loc='center', y=0.975)
    ymax = 0
    for (i, df) in enumerate(pmf_list):
        assert len(df) == 1
        lbl = llist[i]
        col = clist[i]
        df = df[0]
        ymax = max(df.max(), ymax)
        plt.bar(xi+s, df, w, label=lbl, color=col, align='edge')
        s += w

    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.02)
    plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck], y=0.045)
    plt.xlabel('deleterious fitness effect', labelpad=3)
    plt.legend(loc='upper left', ncol=ncol)

    # plot udel
    ax2 = plt.subplot(165)
    format_panels(ax2)
    plt.title('(ii)', loc='center', y=0.975)
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    s = -w / 2.0
    li = 0
    for (i, df) in enumerate(pmf_list):
        for d in df:
            udel = sum(d)
            plt.bar(0+s, udel, w, color=clist[li])
            li += 1
            s += w
    plt.xticks([])
    plt.yticks(x=0.15)

    # plot pi0
    ax3 = plt.subplot(166)
    format_panels(ax3)
    plt.title('(iii)', loc='center', y=0.975)
    plt.ylabel('average reduction\nin heterozygosity', labelpad=3)
    plt.yticks(x=0.15)
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0 + s, p0, w, color=clist[i])
        s += w
    plt.xticks([])

    f_save = final_dir + '/sfigs/cadd_vs_filter_CG_params.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


# PLOT R^2 AND OBS/PRED FOR CADD VS. CADD FILTER CG
def plot_rsq_and_sortpred(flist, clist, llist):
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.995   , right=0.995, hspace=0.1, bottom=0.11,
                        left=0.07, wspace=0.25)
    nbins = 100
    axmin, axmax = 0.52, 1.12
    tx1, tx2 = 0.075, 0.59
    ty1 = 0.943

    # PLOT SORTED RESULT
    ax1 = plt.subplot(121)
    format_panels(ax1)
    for i in xrange(len(flist)):
        f_sort = final_dir + '/{}/sort_gc_cm_cn_n{}.txt'.format(flist[i], nbins)
        div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
        rst = cst_from_fldr(flist[i])
        pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
        pi /= pi0
        pred /= pi0
        # plot predicted vs. observed
        plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
                 markeredgecolor=clist[i], markeredgewidth=0.9, lw=0,
                 alpha=0.75, label=llist[i])
        # plot y=x line
        plt.plot([axmin, axmax], [axmin, axmax], color='darkslategray',
                 ls='--', alpha=0.65)

    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=2)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=2)
    xtick = np.arange(0.5, 1.2, 0.1)
    plt.xticks(xtick, y=0.02)
    plt.ylim(axmin, axmax)
    plt.xlim(axmin, 1.02)
    plt.legend(loc='lower right')
    plt.text(0.75, 0.81, r'$y=x$', rotation=42, ha='center', va='center',
             color='darkslategray', alpha=0.65)
    plt.text(tx1, ty1, 'B', transform=plt.gcf().transFigure)

    # PLOT R^2 RESULT
    ax2 = plt.subplot(122)
    format_panels(ax2)
    for i in xrange(len(flist)):
        f_rsq = final_dir + '/{}/rsq.log'.format(flist[i])
        w, r = np.loadtxt(f_rsq)[:16].T
        plt.plot(np.log10(w), r, marker='o', ms=5, markerfacecolor='None',
                 markeredgecolor=clist[i], markeredgewidth=0.9, lw=0,
                 alpha=0.75)
    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.1f}$' % x if not x % 1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.02)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $\mathrm{(R^2)}$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1,
               labelspacing=0.25)
    plt.text(tx2, ty1, 'C', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/CADD_filter_CG_sort_rsq.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


flist = ['cadd94_gmask_mnb_378', 'cadd94_gmask_filter_CG']
clist = ['darkorange', 'purple']
llist = ['CADD', 'CADD (filter CG)']
compare_params(flist, llist, clist)
# plot_rsq_and_sortpred(flist, clist, llist)


#%% LOAD ARCHAIC INTROGRESSION LEVELS SORTED BY B
def get_archaic(fldr, pop, num):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_arch = fdir + 'predsort_archaic_{}_n{}.txt'.format(pop, num)
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pred, ai = np.loadtxt(f_arch).T
    pred /= pi0

    return pred, ai


def archaic_introgression_plot(fldr, pop, ttl, num, letter):
    pred, ai = get_archaic(fldr, pop, num)
    fig_dir = root_dir + '/result/final_files/sfigs/'
    f_save = fig_dir + 'fig_S44.{}_archaic_{}_n{}.png'.format(fldr, ttl, num)
    plt.figure(figsize=(3.25,3.25))
    plt.subplots_adjust(top=0.93, right=0.99, hspace=0.1, bottom=0.14,
                        left=0.25)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.title(ttl, y=0.98)
    plt.plot(pred, ai, marker='o', ms=5, lw=0, alpha=0.75, color='forestgreen')
    # plt.plot(pred, ai, alpha=0.75, color='forestgreen')

    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel('probability of archaic ancestry', labelpad=3)
    plt.xlim(0.55, 1.02)
    plt.text(0.26, 0.88, letter, transform=plt.gcf().transFigure)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_mnb_378'
n = 100
pops = ['YRI', 'MSL', 'CEU', 'CHBS']
titles = ['African-1', 'African-2', 'European', 'East-Asian']
letters = ['A', 'B', 'A', 'B']
for i in [2, 3]:
    archaic_introgression_plot(fldr, pops[i], titles[i], n, letters[i])
#%% CHECK NEW SORTED CONS DATA
sdir = root_dir + '/compress/sort_cons'
fload = sdir + '/chr22.cadd94_gmask_mnb_378.sort_cons.npz'
a = np.load(fload)['arr']
#%%
bsize = len(a) / 1000
aav = []
for i in xrange(0, len(a), bsize):
    aav.append(np.average(a[i:i+bsize], axis=0))
#%%
aav = np.array(aav)
si = np.argsort(aav[:,0])
#%%
plt.plot(aav[si][:,0], aav[si][:,2])
plt.show()
#%%
plt.plot(a[:,1])
plt.show()
#%%