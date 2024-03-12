__author__ = 'davidmurphy'

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from figures.common_functions import format_panels
from data_processing.functions import relative_error
from classes.runstruct import root_dir, cst_from_fldr
from phast.branch_lengths import mean_matrix, sum_branches

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'


#%% ALIGNMENT COVERAGE TO HUMAN
def alignment_coverage():
    bigfont = 24
    rdir = root_dir + '/result'
    fname = rdir + '/100way_species_coverage/autosomes.100way.cov.txt'
    auto_len = 2881033286.0  # total autosome length

    # get labels and alignment-to-human counts
    with open(fname, 'r') as f:
        a = [l.split() for l in f.readlines()]
        s, n = [], []
        for l in a:
            s.append(' '.join(l[:-1]))
            n.append(float(l[-1]))
        n = np.array(n)

    plt.figure(figsize=(15, 7.5))
    # plt.figure(figsize=(6.5, 3.25))

    plt.subplots_adjust(right=0.995, left=0.08, bottom=0.22, top=0.98)
    ax = plt.subplot(111)
    format_panels(ax)
    plt.bar(left=range(100), height=n / auto_len, color='seagreen')

    # add lines dividing each alignment depth node
    plt.axvline(x=4.5, ls='--')
    plt.text(3.73, 0.96, 'ape', ha='center', va='top', rotation=90, fontsize=18)
    plt.axvline(x=8.5, ls='--')
    plt.text(7.73, 0.96, 'primate', ha='center', va='top', rotation=90, fontsize=18)
    plt.axvline(x=12.5, ls='--')
    plt.text(11.73, 0.96, 'prosimian', ha='center', va='top', rotation=90, fontsize=18)
    plt.axvline(x=25.5, ls='--')
    plt.text(24.73, 0.96, 'supra-primate', ha='center', va='top', rotation=90, fontsize=18)
    plt.axvline(x=50.5, ls='--')
    plt.text(49.73, 0.96, 'laurasiatheria', ha='center', va='top', rotation=90, fontsize=18)
    plt.axvline(x=61.5, ls='--')
    plt.text(60.73, 0.96, 'mammal', ha='center', va='top', rotation=90, fontsize=18)
    plt.axvline(x=99.5, ls='--')
    plt.text(98.73, 0.96, 'fish', ha='center', va='top', rotation=90, fontsize=18)

    # add labels to plot
    plt.xticks(range(100), s, rotation=90, fontsize=10)
    plt.yticks(np.arange(0, 1.01, 0.05), fontsize=bigfont)
    plt.ylabel('fraction aligned to hg19', fontsize=bigfont)
    plt.xlim(-1, 100)
    plt.ylim(0, 0.965)
    plt.text(0.01, 0.95, 'b', fontsize=bigfont, transform=plt.gcf().transFigure,
             fontweight='bold')

    # save figure
    f_save = rdir + '/final_files/sfigs/fig_S10A.autosomes.100way.coverage.png'
    plt.savefig(f_save, dpi=512)
    # plt.show()


alignment_coverage()


# CUMULATIVE PHASTCONS SCORES
def cumulative_phast_scores():
    bigfont = 24
    # set paths to unfiltered file & filtered files
    pth = '/Users/davidmurphy/GoogleDrive/linked_selection/data/cons'
    cols = 'salmon indianred firebrick darkturquoise teal goldenrod ' \
           'royalblue'.split()
    spec = 'ape,primate,prosimian,euarchontoglires,laurasiatheria,' \
           'mammal,fish'.split(',')
    nspc = [4,8,12,25,50,61,99]

    plt.figure(figsize=(15, 7.5))
    plt.subplots_adjust(right=0.995, left=0.08, top=0.995, bottom=0.15)
    ax = plt.subplot(111)
    format_panels(ax)

    # plt.title('cumulative distribution of phastcons scores')
    for i, cons in enumerate(spec):
        f_1 = pth + '/{}.scores.dist.txt'.format(cons)

        # load dist data
        p_1, c_1 = np.loadtxt(f_1).T

        # convert counts to cumulative fraction
        fr_1 = np.cumsum(c_1).astype(float)
        fr_1 /= fr_1[-1]

        # plot cumulant
        lbl = '{} ({})'.format(cons, nspc[i])
        if cons == 'euarchontoglires':
            lbl = 'supra-primate ({})'.format(nspc[i])

        plt.step(p_1, fr_1, color=cols[i], label=lbl, alpha=0.85)

    plt.xlabel('phastcons score', fontsize=bigfont)
    xt = np.arange(0,1.01,0.05)
    plt.xticks(xt, ['{:.2f}'.format(x) for x in xt], rotation=90, fontsize=bigfont)
    plt.ylabel('cumulative fraction', fontsize=bigfont)
    plt.yticks(np.arange(0, 1.01, 0.05), fontsize=bigfont)
    plt.ylim(-0.025, 1.025)
    plt.xlim(-0.01, 1.01)
    plt.legend(prop=dict(size=20))
    plt.text(0.01, 0.96, 'a', fontsize=bigfont, transform=plt.gcf().transFigure,
             fontweight='bold')
    f_save = root_dir + '/result/final_files/sfigs/fig_S10B.combined.scoredist.png'
    plt.savefig(f_save, dpi=256)
    plt.close()
    # plt.show()


cumulative_phast_scores()


#%% COMPARE R^2 AT DIFFERENT WINDOWS FOR DIFFERENT CM MASKED
def rsq_3_size_gmask():
    # gmask file suffixes
    suffs = '000 002 004 006 008 01 02 04 06 08 10'.split()
    # cM equivalents
    cm = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1]
    # files list
    fstring = 'cadd94_gmask_v1.6_without_bstat_cutoff_gmask_{}'
    flist = [fstring.format(s) for s in cm]
    # create string format for rsq files
    fstr = final_dir + '/{}/rsq.log'

    # keep just 3 window sizes: 62.5K, 250K, 1M
    keep_wins = [6.25e4, 2.5e5, 1e6]
    wlab = ['62.5 kb', '250 kb', '1 Mb']
    rsq_list = []
    for fl in flist:
        # load rsq and window size
        win, rsq = np.loadtxt(fstr.format(fl)).T
        # keep rsq at the subset of windows needed
        idx = np.in1d(win, keep_wins)
        rsq_list.append(rsq[idx])

    # convert rsq values to array of columns for each size
    r = np.array(rsq_list)

    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.105, bottom=0.21, wspace=0.3, right=0.99,
                        top=0.88)

    # plot subplot for each window size
    for pidx in range(1, 4):
        ax = plt.subplot(1, 3, pidx)
        format_panels(ax)
        # get column for current window size
        ri = r[:,pidx-1]
        print ri
        # rmin = int(1000*min(ri)) * 1e-3
        # rmax = int(1000*max(ri)+1) * 1e-3
        rmin = round(min(ri)-0.0005, 3)
        rmax = round(max(ri)+0.0005, 3)
        # rmin, rmax = round(min(ri), 3), round(max(ri), 3)
        r_interval = (rmax - rmin) / 5.0
        ytck = np.arange(rmin, rmax+r_interval, r_interval)
        # ytck = [round(rv, 3) for rv in ]

        plt.title(wlab[pidx-1] + ' windows', y=0.96)
        plt.plot(cm, ri, lw=0.8, marker='o', ms=3, color='k',
                 alpha=0.75)
        # plt.plot(0.1, ri[5], lw=0, marker='o', ms=3, color='r')
        # plt.axvline(0.1, color='red', ls='--', lw=0.5, alpha=0.75,
        #             label='0.1 cM')
        plt.xlabel('cM removed from edge', labelpad=1)
        # xtck = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2]
        xtck = np.arange(0, 1.01, 0.1)
        xtck_show = [0.1, 0.4, 1]
        plt.xticks(xtck, [x if x in xtck_show else '' for x in xtck])
        # plt.xticks(xi, [l[:9] for l in llist], rotation=rotation, y=0.06,
        #            ha='center')
        plt.yticks(ytck, x=0.05)
        plt.ylim(rmin - 0.5*r_interval, rmax + 0.5*r_interval)
        if pidx == 1:
            ylab = r'variance explained ($R^2$)'
            plt.ylabel(ylab, labelpad=1)
            # plt.legend(frameon=1, framealpha=0.75, facecolor='white',
            #            loc='upper right')
    # figure letter
    plt.text(0.01, 0.95, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')
    # plt.text(0.43, 0.81, 'B) 250Kb windows', transform=plt.gcf().transFigure)
    # plt.text(0.75, 0.81, 'C) 1Mb windows', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/fig_S12.gmask.3range.rsq.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


#   cadd93_gmask_002
# cadd93_gmask_

rsq_3_size_gmask()


#%% COMPARE GMASKING INFERRED PARAMS
def compare_gmask_params():
    # gmask file suffixes
    suffs = '000 002 004 006 008 01 02 04 06 08 10'.split()
    # cM equivalents
    cm = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1]
    # files list
    # cadd94_gmask_v1.6_without_bstat_cutoff_gmask_0.4
    fstring = 'cadd94_gmask_v1.6_without_bstat_cutoff_gmask_{}'
    flist = [fstring.format(s) for s in cm]
    prm_idx = [0, 2, 4, 5, 7, -2, -1]
    # colors list
    clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue',
             'rosybrown', 'darkturquoise']
    cm_list = [cm[i] for i in prm_idx]
    xtck_show = [0.1, 0.4, 1]
    xtck_str = '0.1 0.4 1.0'.split()
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
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.07, wspace=0.55, right=0.985, bottom=0.17,
                        top=0.91)
    alpha_val = 0.4
    shaded_col = 'white'

    # DFE plot
    ax1 = plt.subplot(1, 4, (1, 2))
    format_panels(ax1)

    plt.title('b', loc='left', y=0.975, fontweight='bold')
    plt.title('(i)', loc='center', y=0.975, fontweight='bold')
    ymax = 0
    # plt.text(0.08, 0.83, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.08, 0.83, '(shaded = phastCons)', transform=plt.gcf().transFigure)
    # plt.text(0.08, 0.73, '(solid = CADD)', transform=plt.gcf().transFigure)
    # hatches = ['.', '/', 'x']
    j = 0
    for i in prm_idx:
    # for (i, df) in enumerate(pmf_list):
        df = pmf_list[i]
        assert len(df) == 1
        lbl = str(cm[i])+'cM'
        col = clist[j]
        df = df[0]
        ymax = max(df.max(), ymax)
        plt.bar(xi+s, df, w, label=lbl, color=col, align='edge')
        s += w
        j += 1

    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    plt.ylim(0, 1.5*ymax)
    plt.yticks(x=0.02)
    plt.xticks(xi, xtck, y=0.05)
    plt.xlabel(r'deleterious fitness effect ($\mathrm{log_{10}}$ scale)',
               labelpad=1)
    plt.legend(loc='lower left', ncol=ncol)

    # plot udel
    ax2 = plt.subplot(143)
    format_panels(ax2)
    plt.title('(ii)', loc='center', y=0.975, fontweight='bold')
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    ulist = []
    for i in prm_idx:
        df = pmf_list[i]
        for d in df:
            udel = sum(d)
            ulist.append(udel)
    plt.plot(cm_list, ulist, marker='o', ms=3, color='k', lw=0.8, alpha=0.75)
    xtck = np.arange(0, 1.01, 0.1)
    xtck_show = [0.1, 0.4, 1]
    plt.xticks(xtck, [x if x in xtck_show else '' for x in xtck])
    # plt.xticks(xtck_show, xtck_str, y=0.05)
    plt.xlabel('cM removed', labelpad=1)
    plt.yticks(x=0.08)

    # plot pi0
    ax3 = plt.subplot(144)
    format_panels(ax3)
    plt.title('(iii)', loc='center', y=0.975, fontweight='bold')
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=1)
    plt.plot(cm_list, [pi0_list[i] for i in prm_idx], marker='o', ms=3,
             color='k', lw=0.8, alpha=0.75)
    plt.yticks(x=0.08)
    plt.xticks(xtck, [x if x in xtck_show else '' for x in xtck])

    plt.xlabel('cM removed', labelpad=1)

    # plt.text(0.01, 0.95, 'b', transform=plt.gcf().transFigure,
    #          fontweight='bold')

    f_save = final_dir + '/sfigs/gmask.params.png'
    plt.savefig(f_save, dpi=512)


compare_gmask_params()


#%% NEUTRAL CHOICES (NEUT ALIGNMENT AND NEUT CUTOFF SCORE)
def neutral_choices_1():
    """rsq for different neutral choices"""
    spec = 'ape primate prosimian euarchontoglires ' \
           'laurasiatheria mammal fish'.split()
    # val = range(7)
    nvals = [0, 1, 2, 3, 4, 5, 10, 20, 30]

    # set chr1 file name
    s_file = final_dir + '/species_rsq.txt'
    # p_file = final_dir + '/percent_rsq.txt'
    w_file = final_dir + '/window_rsq.txt'
    v_file = final_dir + '/values_rsq.txt'

    # load R^2 data for different species alignments
    s_fish = {}
    s_cad = {}
    with open(s_file, 'r') as f:
        for line in f:
            fldr, r2 = line.split()
            r2 = float(r2)
            sp = fldr.split('_')[-1]
            if fldr.startswith('fish'):
                s_fish[sp] = r2
            else:
                assert fldr.startswith('cadd')
                s_cad[sp] = r2

    # load R^2 data for different neutral sub window sizes
    w_fish = {}
    w_cad = {}
    with open(w_file, 'r') as f:
        for line in f:
            fldr, r2 = line.split()
            r2 = float(r2)
            wn = str(int(fldr.split('_')[-1][3:]) / 1000)
            if fldr.startswith('fish'):
                w_fish[wn] = r2
            else:
                assert fldr.startswith('cadd')
                w_cad[wn] = r2

    # load R^2 data for different phastcon value cutoffs
    v_fish = {}
    v_cad = {}
    with open(v_file, 'r') as f:
        for line in f:
            fldr, r2 = line.split()
            r2 = float(r2)
            # vl = int(fldr[-1])
            vl = int(fldr.split('_')[-1][4:])
            # if fldr.startswith('ape'):
            if fldr.startswith('fish'):
                v_fish[vl] = r2
            else:
                assert fldr.startswith('cadd')
                v_cad[vl] = r2

    # create new figure
    plt.figure(figsize=(5.4, 2.16))
    plt.subplots_adjust(left=0.13, bottom=0.2, right=0.97, top=0.98,
                        wspace=0.38)

    # plot species R^2
    ax1 = plt.subplot(121)
    format_panels(ax1)
    xi = range(len(spec))
    plt.plot(xi, [s_fish[sp] for sp in spec], color='darkorange',
             marker='o', lw=0, ms=7, alpha=0.8,label='phastCons')
    plt.plot(xi, [s_cad[sp] for sp in spec], color='dodgerblue',
             marker='o', lw=0, ms=7, alpha=0.8, label='CADD')

    plt.xticks(xi, [sp[:4] for sp in spec], rotation=45, y=0.03)
    plt.xlabel('alignment', labelpad=2)
    plt.yticks(x=0.04)
    plt.ylim(0.43, 0.62)
    plt.ylabel('variance explained in\n1Mb windows ' + r'($R^2$)', labelpad=2)
    # plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1)

    # plot value R^2
    ax2 = plt.subplot(122)
    format_panels(ax2)

    xi = range(len(nvals))
    plt.plot(xi, [v_fish[vl] for vl in nvals], color='darkorange',
             marker='o', lw=0, ms=7, alpha=0.8,label='phastCons')
    plt.plot(xi, [v_cad[vl] for vl in nvals], color='dodgerblue',
             marker='o', lw=0, ms=7, alpha=0.8, label='CADD')

    plt.xticks(xi, [vl for vl in nvals], rotation=45, y=0.03)
    plt.xlabel(r'cutoff phastCons score $(\times 1000)$', labelpad=2)
    plt.ylabel('variance explained in\n1Mb windows ' + r'($R^2$)', labelpad=2)
    # plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks(x=0.04)
    plt.ylim(0.43, 0.62)

    # plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    # plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1)
    plt.text(0.137, 0.91, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.625, 0.91, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/supp_3_rsq_neutral_choices.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


neutral_choices_1()


#%% BRANCH LENGTH VS. FRACTION ALIGNED IN 8 PRIMATES
def branchlength_vs_fractionaligned():
    fmt = root_dir + '/data/phast/primate_neutral/{}.neutral.exptotsub'
    astats = root_dir + '/data/phast/align_stats_combined.txt'
    fsave = root_dir + '/result/final_files/sfigs/fig_S13.branch_lengths.png'
    clist = 'darkturquoise darkorange purple gold fuchsia limegreen blue firebrick'.split()

    mm = mean_matrix(fmt)
    lpct, lcdist, lsdist, lnum = [], [], [], []

    with open(astats, 'r') as f:
        for line in f:
            sp, _, pct = line.split()
            slist = sp.split('_')[1:]
            cdist, sdist = sum_branches(slist, mm)
            lpct.append(float(pct))
            lnum.append(len(slist))
            lcdist.append(cdist)
            lsdist.append(sdist)

    lpct, lcdist, lsdist, lnum = [np.array(a) for a in lpct, lcdist, lsdist, lnum]

    # max n*p
    max_np = max(l * p for (p, l) in zip(lpct, lsdist))
    # print the length associated with max n*p
    print 'max length',[l for (p, l) in zip(lpct, lsdist) if l*p == max_np][0]

    pdict = defaultdict(list)
    ldict = defaultdict(list)
    for (p, l, n) in zip(lpct, lsdist, lnum):
        pdict[n].append(p)
        ldict[n].append(l)

    plt.figure(figsize=(6.5, 3.25))
    ax1 = plt.subplot(121)
    format_panels(ax1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.995, top=0.98,
                        wspace=0.38)
    for n in sorted(pdict.keys()):
        plt.scatter(pdict[n], ldict[n], color=clist[n-1], label=str(n+1),
                    alpha=0.8)
    # line of equal n*p for max np
    xi = np.arange(0.8, 0.99, 0.001)
    fx = lambda x: max_np / x
    plt.plot(xi, fx(xi), color='gray', ls='--')
    plt.text(0.94, 0.15, r'$max(n \cdot T)$', color='grey', ha='center')

    plt.xlabel('fraction of aligned sites', labelpad=2)
    plt.ylabel('total branch length (substitutions per site)', labelpad=2)
    plt.xlim(0.8, 0.99)
    plt.xticks(np.arange(0.8, 1.0, 0.04), y=0.02)
    plt.yticks(x=0.03)
    plt.ylim(0.015, 0.185)
    plt.legend(title='# species', ncol=2, loc='lower left', borderpad=0.2,
               frameon=1, framealpha=0.75, facecolor='white',
               columnspacing=0.1, labelspacing=0.25)

    # plot window size R^2
    ax2 = plt.subplot(122)
    format_panels(ax2)
    # load R^2 data for different neutral sub window sizes
    w_fish = {}
    w_cad = {}
    w_file = final_dir + '/window_rsq.txt'
    win = '4 5 6 7 8 9 10 11'.split()
    # scale_factor = 0.59746/0.59575
    with open(w_file, 'r') as f:
        for line in f:
            fldr, r2 = line.split()
            r2 = float(r2)
            wn = str(int(fldr.split('_')[-1][3:]) / 1000)
            # if fldr.startswith('ape'):
            if fldr.startswith('fish'):
                w_fish[wn] = r2
            else:
                assert fldr.startswith('cadd')
                w_cad[wn] = r2

    xi = range(len(win))
    plt.plot(xi, [w_fish[wn] for wn in win], color='darkorange',
             marker='o', lw=0, ms=7, alpha=0.8, label='phastCons')
    plt.plot(xi, [w_cad[wn] for wn in win], color='dodgerblue',
             marker='o', lw=0, ms=7, alpha=0.8, label='CADD')

    plt.xticks(xi, [wn for wn in win], y=0.02)
    plt.xlabel('window size (in kb)', labelpad=2)
    plt.ylim(0.57, 0.6)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks(x=0.03)
    plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
               borderaxespad=0.3, ncol=1, frameon=1,
               framealpha=0.75, facecolor='white', columnspacing=0.1)

    plt.text(0.105, 0.925, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.625, 0.925, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.savefig(fsave, dpi=512)
    plt.close()


branchlength_vs_fractionaligned()


#%% STATS ON WINDOW SIZES (REL ERR, STD, MEAN)
def window_size_mean_and_std(win_size):
    """get mean & std for PHYSICAL window sizes from constant count windows"""
    ddir = root_dir + '/data/phast/aligned_8_win_{}'.format(win_size)

    # get length of each window for each chrom
    windows = []
    n, k = [], []
    st, en = [], []
    for c in range(1, 23):
        f_in = ddir + '/chr{}.subcount.filled.txt'.format(c)
        start, end, win, subs = np.loadtxt(f_in, usecols=(0,1,2,3)).T
        size = end-start
        st.append(start)
        en.append(end)
        windows.append(size)
        n.append(win)
        k.append(subs)

    windows = np.concatenate(windows)
    n, k = np.concatenate(n), np.concatenate(k)
    # st, en = np.concatenate(st), np.concatenate(en)

    print 'mean mu = {}'.format(np.mean(k/n))
    # calculate mean and std of windows
    wmean, wstd = np.mean(windows), np.std(windows)

    # calculate relative error in mutation rate estimate
    relerr = relative_error(n, k)

    # get the standard error of relative errors
    serelerr = np.std(relerr) / np.sqrt(len(relerr))

    return wmean, wstd, np.mean(relerr), serelerr, windows


wmean, wstd, werr, sererr, wins = window_size_mean_and_std(6000)
#%%