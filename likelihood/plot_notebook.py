__author__ = 'davidmurphy'


import os, re
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import t as ttest
from collections import defaultdict
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir, \
    human_autosomes, RunStruct
from classes.phylotree import load_tree, parse_tree
from phast.branch_lengths import mean_matrix
import matplotlib_venn as venn
# dfe = default_fitness_effects
from scipy.stats import pearsonr, spearmanr


def compare_runs(folder_list, label_list, colors, markers, save_name):
    # create lists for data to plot
    dfe_list = []
    pmf_list = []
    pi0_list = []
    clh_list = []
    rsq_list = []

    for fldr in folder_list:
        # find composite file and rsq file in folder (should be 1 of each)
        fdir = root_dir + '/result/final_files/{}/'.format(fldr)
        flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
        rlist = [f for f in os.listdir(fdir) if f.endswith('.log')]
        assert len(flist) == len(rlist) == 1

        # create paths to rsq and final init file
        f_init = fdir + flist[0]
        r_log = fdir + rlist[0]

        # get the rsq values
        rsq = np.loadtxt(r_log)
        rsq_list.append(rsq)

        # initalize ChromStruct with init file for final params
        cst = ChromStruct('chr1', init=f_init)

        # get dfe and pi0 from composite run
        # mpmf = cst.uvec[0] * 1e8
        mpmf = [u*1e8 for u in cst.uvec]
        mpi0 = cst.params[-1] / cst.fixed.tau_init
        dfe_list.append(cst.bdfe[0])
        pmf_list.append(mpmf)
        pi0_list.append(mpi0)

        # get best CLH
        clh_list.append(cst.stat.best_lh)

    # x = np.arange(4)
    x = np.arange(6)

    n = len(folder_list)
    w = 0.8 / n
    s = -w * n / 2.0

    # create the plot
    plt.figure(figsize=(12, 4))
    plt.suptitle(save_name.replace('_', ' '))
    plt.subplots_adjust(left=0.05, wspace=1.5, right=0.95)

    # plot DFE
    if any(len(pm) > 1 for pm in pmf_list):
        plt_1 = plt.subplot(2, 11, (1, 4))
        plt_2 = plt.subplot(2, 11, (12, 15))
        ymax = 0
        for (i, df) in enumerate(pmf_list):
            lbl = label_list[i]
            if len(df) > 1:
                ymax = max(max(d.max() for d in df), ymax)
                d = df[0]
                # fix for DFE with truncation
                if len(d) < 6:
                    nz = 6 - len(d)
                    zpad = [0] * nz
                    d = np.concatenate((zpad, d))
                plt_1.bar(x + s, d, w, label=lbl, color=colors[i],
                          align='edge')
                d = df[1]
                # fix for DFE with truncation
                if len(d) < 6:
                    nz = 6 - len(d)
                    zpad = [0] * nz
                    d = np.concatenate((zpad, d))
                plt_2.bar(x + s, d, w, label=lbl, color=colors[i], align='edge')
            else:
                plt_1.bar(x + s, df[0], w, label=lbl, color=colors[i],
                          align='edge')
            s += w
        plt_1.set_title('distribution of fitness effects')
        plt_1.set_ylabel(r'$\mu_{del} \cdot 10^8$')
        plt_2.set_ylabel(r'$\mu_{del} \cdot 10^8$')
        plt_1.set_xticks([])
        plt.ylim(0, 1.5 * ymax)
        xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]

        # xtck = [-3.5, -3.0, -2.5, -2.0]
        plt.xticks(x, [r'$10^{%.1f}$' % x for x in xtck])
        plt.xlabel('deleterious fitness effect')
        plt_1.legend(prop=dict(size=9), loc='upper left')
        plt_2.legend(['conserved nonexonic'], prop=dict(size=9),
                     loc='upper left')
    else:
        plt.subplot(1, 11, (1, 4))
        plt.title('distribution of fitness effects')
        ymax = 0
        for (i, df) in enumerate(pmf_list):
            lbl = label_list[i]
            df = df[0]
            # fix for DFE with truncation
            if len(df) < 6:
                nz = 6 - len(df)
                zpad = [0] * nz
                df = np.concatenate((zpad, df))
            ymax = max(df.max(), ymax)
            plt.bar(x + s, df, w, label=lbl, color=colors[i], align='edge')
            s += w
        plt.ylabel(r'$\mu_{del} \cdot 10^8$')
        plt.ylim(0, 1.5*ymax)
        # xtck = [-3.5, -3.0, -2.5, -2.0]
        xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]

        plt.xticks(x, [r'$10^{%.1f}$' % x for x in xtck])
        plt.xlabel('deleterious fitness effect')
        plt.legend(prop=dict(size=9), loc='upper left')

    # # plot dfe higher
    # x = np.arange(3)
    # n = 2
    # w = 0.8 / n
    # s = -w * n / 2.0
    # plt.subplot(2, 11, (1, 4))
    # plt.title('distribution of fitness effects')
    # for (i, df) in enumerate(pmf_list[:2]):
    #     lbl = label_list[i]
    #     plt.bar(x + s, df, w, label=lbl, color=colors[i], align='edge')
    #     s += w
    # plt.ylabel(r'$\mu_{del} \cdot 10^8$')
    # xtck = [-3.0, -2.5, -2.0]
    # plt.xticks(x, [r'$10^{%.1f}$' % x for x in xtck])
    # # plt.xlabel('deleterious fitness effect')
    # plt.legend(prop=dict(size=9))
    #
    # # plot dfe lower
    # x = np.arange(6)
    # n = 2
    # w = 0.8 / n
    # s = -w * n / 2.0
    # plt.subplot(2, 11, (12, 15))
    # # plt.title('distribution of fitness effects')
    # for (i, df) in enumerate(pmf_list[2:]):
    #     lbl = label_list[i+2]
    #     plt.bar(x + s, df, w, label=lbl, color=colors[i+2], align='edge')
    #     s += w
    # plt.ylabel(r'$\mu_{del} \cdot 10^8$')
    # xtck = [-3.0, -2.8, -2.6, -2.4, -2.2, -2.0]
    # plt.xticks(x, [r'$10^{%.1f}$' % x for x in xtck])
    # plt.xlabel('deleterious fitness effect')
    # plt.legend(prop=dict(size=9))

    # plot udel
    sp3 = plt.subplot(1,11,5)
    sp3.yaxis.tick_right()
    sp3.yaxis.set_label_position('right')
    plt.xlabel(r'$\mu_{del} \cdot 10^8$')
    s = -w / 2.0
    for (i, df) in enumerate(pmf_list):
        if len(df) > 1:
            # udel = 0.58*(sum(df[0]) + 0.42*sum(df[1]))
            df = [1.5305702927883065, 0.8530486891064619]
            udel = 0.12 * df[0] + 0.88 * df[1]
        else:
            udel = sum(df[0])
        plt.bar(0+s, udel, w, color=colors[i])
        s += w
    plt.xticks([])

    # plot pi0
    sp4 = plt.subplot(1,11,6)
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    plt.xlabel(r'$\frac{\pi}{\pi_0}$')
    s = -w / 2.0
    for (i, p0) in enumerate(pi0_list):
        plt.bar(0+s, p0, w, color=colors[i])
        s += w
    plt.xticks([])
    # plt.yticks([])

    # plot CLH
    sp4 = plt.subplot(1,11,7)
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    sp4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('-CLLH')
    s = -w * n / 2.0
    clh_list = [cl/1e5 for cl in clh_list]
    for (i, cl) in enumerate(clh_list):
        print '{:<20} -CLLH {}'.format(folder_list[i], cl)
        plt.bar(0+s, cl, w, color=colors[i])
        s += w
    plt.ylim(min(clh_list), max(clh_list))
    # plt.ylim(0.999*min(clh_list), 1.001*max(clh_list))
    plt.xticks([])

    sp4 = plt.subplot(1,11,(8,11))
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    plt.title('variance explained')
    plt.ylabel(r'$R^2$')
    # get window sizes from the first file
    ws = np.log10(rsq_list[0][:,0])
    wmsk = (ws < 6.4)
    # s = -w * n / 2.0
    for (i, rs) in enumerate(rsq_list):
        rsq = rs[:,1]
        plt.plot(ws[wmsk], rsq[wmsk], color=colors[i], marker=markers[i],
                 linewidth=0, alpha=0.75, label=label_list[i])
        # plt.bar(x+s, rsq, w, color=colors[i], align='edge')
        # s += w
    plt.xlabel('window size')
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4,5,6]])
    plt.legend(prop=dict(size=9), loc='upper left')

    # save the plot
    save_lbl = 'compare_cons_'+'_'.join(save_name.split())
    fsave = root_dir + '/result/final_files/{}.png'.format(save_lbl)
    plt.savefig(fsave, dpi=256)
    plt.close()
    # plt.show()


def plot_observed_vs_expected_cllh(folders, labels, colors, figname):
    # load data

    arrs = []
    for fldr in folders:
        # load list of cllh gap array
        fa = root_dir + '/result/final_files/{}/cllh_gap_100.txt'.format(fldr)
        arrs.append(np.loadtxt(fa, dtype=float))

    # plot indexed results
    plt.figure(figsize=(10, 5))
    for (i, arr) in enumerate(arrs):
        bi, ob, ex = arr.T
        plt.plot(ob/ex, label=labels[i], color=colors[i], alpha=0.9, lw=0.8)
    plt.ylabel('observed/expected CLLH')
    plt.xlabel('data bin')
    plt.legend()
    fsave = root_dir + '/result/final_files/{}_idx.png'.format(figname)
    plt.savefig(fsave, dpi=256)
    plt.close()

    # plot indexed results
    plt.figure(figsize=(10, 5))
    for (i, arr) in enumerate(arrs):
        bi, ob, ex = arr.T
        plt.plot(bi, ob/ex, label=labels[i], color=colors[i], alpha=0.9, lw=0.8)
    plt.axvline(x=0.65, ls='--', alpha=0.25, color='k')
    plt.ylabel('observed/expected CLLH')
    plt.xlabel('B value bin')
    plt.legend()
    fsave = root_dir + '/result/final_files/{}_bx.png'.format(figname)
    plt.savefig(fsave, dpi=256)
    plt.close()


def plot_obex_regression(f_obex):
    # load obex array elements
    bsx, ehom, ehet, ohom, ohet = np.load(f_obex).T

    # mask out inadvertant 0 site rows
    ms = (ehom > 0)
    assert np.all(ms == (ohom > 0))
    # bsx = bsx[ms]
    ehom, ehet, ohom, ohet = ehom[ms], ehet[ms], ohom[ms], ohet[ms]

    # calculate obs/exp pi
    opi = 1.0 * ohet / (ohom + ohet)
    epi = 1.0 * ehet / (ehom + ehet)

    # calculate regression
    slope, intercept, rvalue, pvalue, stderr = linregress(epi, opi)

    # get endpoints for a line using linear regression line
    xvals = [epi.min(), epi.max()]
    yvals = [(xi * slope) + intercept for xi in xvals]

    # plot scatter plot of obs/exp including regression line
    plt.figure()
    plt.title('expected vs. observed heterozygosity in lowest 5% B bins',
              size=14)
    plt.scatter(x=epi, y=opi, alpha=0.5)
    plt.xlim(epi.min(), epi.max())
    plt.plot(xvals, yvals, color='darkgray')
    plt.xlabel('expected heterozygosity', size=12)
    plt.xticks(size=12)
    plt.text(0.003, 0.3, r'$y = {:.3f}x + {:.6f}$'.format(slope, intercept))
    plt.ylabel('observed heterozygosity', size=12)
    # plt.ylabel('observed heterozygosity (log scaled)', size=12)
    # plt.yscale('log')
    plt.yticks(size=12)
    fsave = f_obex.replace('.npy', '.scatter.png')
    # fsave = f_obex.replace('.npy', '.scatter.log.png')
    plt.savefig(fsave, dpi=256)
    plt.close()


def individual_rsq_plot(flist, llist, sname):
    fstr = root_dir + '/result/final_files/{}/rsq.log'
    win, rs = np.loadtxt(fstr.format(flist[0]))[:16].T
    rsq = [rs]
    for fl in flist[1:]:
        rsq.append(np.loadtxt(fstr.format(fl))[:16,1])
    rsq = np.column_stack(rsq)
    kb = np.arange(len(llist))
    fmt = root_dir + '/result/final_files/{}.{}.rsq.png'
    for r, w in zip(rsq, win):
        plt.figure(figsize=(5,5))
        plt.subplots_adjust(left=0.2, bottom=0.25)
        plt.title('window size = {:.1f}kb'.format(w/1000.0))
        plt.plot(kb, r, lw=0, marker='o', ms=10)
        plt.xlabel('B threshold')
        plt.xticks(kb, [l[:8] for l in llist], rotation=90)
        plt.ylabel(r'$R^2$')
        plt.savefig(fmt.format(sname, w))
        plt.close()


def grouped_rsq():
    clist = 'salmon indianred firebrick darkturquoise teal goldenrod ' \
            'royalblue'.split()
    fstr = root_dir + '/result/final_files/{}/rsq.log'
    spc = 'ape,primate,prosimian,euarchontoglires,laurasiatheria,mammal,fish'
    spc = spc.split(',')
    rdic = defaultdict(dict)
    win = []
    for sc in spc:
        for sn in spc:
            f = fstr.format('{}95_{}35'.format(sc, sn))
            win, r = np.loadtxt(f)[:16].T
            rdic[sc][sn] = r

    kb = np.arange(len(spc))
    fmt = root_dir + '/result/final_files/grouped.{}.rsq.png'
    for (i, w) in enumerate(win):
        plt.figure(figsize=(5,5))
        plt.subplots_adjust(left=0.15)
        plt.title('window size = {:.1f}kb'.format(w/1000.0))
        for (j, sc) in enumerate(spc):
            rsq = [rdic[sc][sn][i] for sn in spc]
            plt.plot(kb, rsq, lw=0, marker='o', ms=10, color=clist[j], label=sc)
        plt.xlabel('conserved alignment')
        plt.xticks(kb, [l[:5] for l in spc])
        plt.ylabel(r'$R^2$')
        plt.legend()
        plt.savefig(fmt.format(w))
        plt.close()


def grouped_rsq_vary_neut():
    clist = 'firebrick goldenrod royalblue'.split()
    fstr = root_dir + '/result/final_files/{}/rsq.log'
    cons_1 = '92 95 98'.split()
    cons_2 = [str(n) for n in range(15, 96, 10) + [100]]
    rdic = defaultdict(dict)
    win = []
    for sc in cons_1:
        for sn in cons_2:
            f = fstr.format('fish{}_euarchontoglires{}'.format(sc, sn))
            win, r = np.loadtxt(f)[:16].T
            rdic[sc][sn] = r

    kb = np.arange(len(cons_2))
    fmt = root_dir + '/result/final_files/vary_neut.{}.rsq.png'
    for (i, w) in enumerate(win):
        plt.figure(figsize=(5,5))
        plt.subplots_adjust(left=0.15)
        plt.title('window size = {:.1f}kb'.format(w/1000.0))
        for (j, sc) in enumerate(cons_1):
            rsq = [rdic[sc][sn][i] for sn in cons_2]
            plt.plot(kb, rsq, lw=0, marker='o', ms=10, color=clist[j], label=sc)
        plt.xlabel('percent neutral cutoff')
        plt.xticks(kb, cons_2)
        plt.ylabel(r'$R^2$')
        plt.legend()
        plt.savefig(fmt.format(w))
        plt.close()


def grouped_rsq_vary_spc_cons_neut(pcon, filt=True):
    clist = 'salmon indianred firebrick darkturquoise teal goldenrod ' \
            'royalblue'.split()
    fstr = root_dir + '/result/final_files/{}/rsq.log'
    # cons_1 = 'ape primate prosimian euarchontoglires laurasiatheria mammal ' \
    #          'fish'.split()
    cons_1 = 'ape euarchontoglires fish'.split()
    cons_2 = [str(n) for n in range(20, 100, 10)]
    if not filt:
        cons_2.append(100)
    rdic = defaultdict(dict)
    win = []
    fmt = '{}{}_euarchontoglires{}'
    if filt:
        fmt += '_filtered'
    for sc in cons_1:
        for sn in cons_2:
            f = fstr.format(fmt.format(sc, pcon, sn))
            win, r = np.loadtxt(f)[:16].T
            rdic[sc][sn] = r

    kb = np.arange(len(cons_2))
    if filt:
        sfmt = root_dir + '/result/final_files/cons{}.filtered.{}.rsq.png'
    else:
        sfmt = root_dir + '/result/final_files/cons{}.{}.rsq.png'
    for (i, w) in enumerate(win):
        plt.figure(figsize=(5,5))
        plt.subplots_adjust(left=0.15)
        plt.title('window size = {:.1f}kb'.format(w/1000.0))
        for (j, sc) in enumerate(cons_1):
            rsq = [rdic[sc][sn][i] for sn in cons_2]
            plt.plot(kb, rsq, lw=0, marker='o', ms=10, color=clist[j],
                     label='{} {}%'.format(sc, pcon))
        plt.xlabel('percent neutral cutoff')
        plt.xticks(kb, cons_2)
        plt.ylabel(r'$R^2$')
        # plt.legend()
        legend = plt.legend(frameon=1)
        frame = legend.get_frame()
        frame.set_color('white')
        plt.savefig(sfmt.format(pcon, w), dpi=256)
        plt.close()


def cumulative_phast_scores(cons):
    # set paths to unfiltered file & filtered files
    pth = '/Users/davidmurphy/GoogleDrive/linked_selection/data/cons/'
    f_1 = pth + '/{}.scores.dist.txt'.format(cons)
    f_2 = pth + '/{}.filtered.scores.dist.txt'.format(cons)

    # load dist data
    p_1, c_1 = np.loadtxt(f_1).T
    # adjust old float scores to 0-1000 scale
    p_1 *= 1000
    p_2, c_2 = np.loadtxt(f_2).T

    # convert counts to cumulative fraction
    fr_1 = np.cumsum(c_1).astype(float)
    fr_1 /= fr_1[-1]
    fr_2 = np.cumsum(c_2).astype(float)
    fr_2 /= fr_2[-1]

    for p in np.arange(0.9, 0.99, 0.01):
        i = np.searchsorted(fr_1, p)
        print '{} {}% = {}'.format(cons, 100*(1-p), p_2[i])

    # plot cumulant
    plt.figure(figsize=(7, 7))
    plt.title(cons)
    plt.step(p_1, fr_1, color='salmon', label='original')
    plt.step(p_2, fr_2, color='teal', label='filtered')
    plt.xlabel('phastcons score')
    plt.xticks(range(0, 1001, 50), rotation=90)
    plt.ylabel('cumulative fraction')
    plt.yticks(np.arange(0, 1.01, 0.05))
    plt.ylim(fr_1.min()-0.025, 1.025)
    plt.legend()
    fsave = pth + '{}.scoredist.png'.format(cons)
    plt.savefig(fsave, dpi=256)
    plt.close()
    # plt.show()


def cumulative_phast_scores_2():
    # set paths to unfiltered file & filtered files
    pth = '/Users/davidmurphy/GoogleDrive/linked_selection/data/cons'
    cols = 'salmon indianred firebrick darkturquoise teal goldenrod ' \
                'royalblue'.split()
    spec = 'ape,primate,prosimian,euarchontoglires,laurasiatheria,' \
    'mammal,fish'.split(',')

    plt.figure(figsize=(7, 7))
    plt.title('cumulative distribution of phastcons scores')
    for i, cons in enumerate(spec):
        f_1 = pth + '/{}.scores.dist.txt'.format(cons)

        # load dist data
        p_1, c_1 = np.loadtxt(f_1).T
        # adjust old float scores to 0-1000 scale
        p_1 *= 1000

        # convert counts to cumulative fraction
        fr_1 = np.cumsum(c_1).astype(float)
        fr_1 /= fr_1[-1]

        # for p in np.arange(0.9, 0.99, 0.01):
        #     i = np.searchsorted(fr_1, p)
        #     print '{} {}% = {}'.format(cons, 100*(1-p), p_2[i])

        # plot cumulant

        plt.step(p_1, fr_1, color=cols[i], label=cons)
        plt.xlabel('phastcons score')
        plt.xticks(range(0, 1001, 50), rotation=90)
        plt.ylabel('cumulative fraction')
        plt.yticks(np.arange(0, 1.01, 0.05))
        plt.ylim(-0.025, 1.025)
        plt.legend()
    fsave = pth + '/combined.scoredist.png'
    plt.savefig(fsave, dpi=256)
    plt.close()


def plot_phylofit_neut_subrates(intervals, sname):
    """plot substitution rates at different intervals of conservation scores"""
    # directories for substitution models and masks
    dmod = root_dir + '/data/phast/euarchontoglires_neutral'
    dmsk = root_dir + '/data/phast/euarchontoglires_masks'
    # file formats for model and mask files
    fmod = dmod + '/{ch}.neut_align.{ii}.{ij}.mod'
    fmsk = dmsk + '/{ch}.euarchontoglires.{ii}.{ij}.nmsk.npz'

    # create dictionaries for interval data, list for interval labels
    intvl_dict = defaultdict(list)
    sum_dict = defaultdict(int)
    labels = []

    # process substitution rates for the first set of intervals
    for (i, j) in intervals:
        # create label for the interval
        lbl = '{}-{}%'.format(i, j)
        labels.append(lbl)
        for ch in human_autosomes:
            # skip chr21, there is no data
            if ch == 'chr21':
                continue

            # set model and mask file names
            mod_file = fmod.format(ch=ch, ii=i, ij=j)
            mask_file = fmsk.format(ch=ch, ii=i, ij=j)

            # get the sum of neutral sites for the chrom
            sites = np.load(mask_file)['neutmask'].sum()
            sum_dict[lbl] += sites

            # check that mod file exists and if not print warning
            if os.path.isfile(mod_file):
                # get the sum of the branch lengths in the tree using regex
                tree_string = load_tree(mod_file)
                subrate = sum(map(float,
                                  re.findall('0\.\d+', tree_string)))

                # add subrate and sites to dict for current label
                intvl_dict[lbl].append((subrate, sites))
            else:
                print 'WARNING: missing model file {}'.format(mod_file)

    # sum site counts and (weighted) average subrates for each label
    counts = []
    rates = []
    for lbl in labels:
        sum_rates = sum(i[0]*i[1] for i in intvl_dict[lbl])
        nsites = sum(i[1] for i in intvl_dict[lbl])
        mean_rate = sum_rates / nsites
        assert nsites == sum_dict[lbl]
        # counts.append(nsites)
        counts.append(sum_dict[lbl])
        rates.append(mean_rate)

    # plot results
    x = range(len(labels))
    plt.figure(figsize=(8, 10))
    plt.subplots_adjust(left=0.15)
    # plot counts on top
    plt.subplot(211)
    plt.bar(x, counts)
    plt.xticks([])
    plt.ylabel('number of neutral sites per bin')
    # plot rates on bottom
    plt.subplot(212)
    plt.bar(x, rates)
    plt.xticks(x, labels)
    plt.ylabel('substitutions per bp in tree')
    figname = root_dir + '/data/phast/{}.png'.format(sname)
    plt.savefig(figname, dpi=256)
    plt.close()
    # plt.show()


def jackknife_estimate(jk_samples):
    """calculate the jackknife estimate from jackknife resampled data"""
    pass


def jackknife_variance(jk_estimates):
    n = len(jk_estimates)
    jk_mean = np.average(jk_estimates, axis=0)
    jk_var = np.sum(np.power(jk_estimates-jk_mean, 2), axis=0) * (n-1) / n

    return jk_var


def jackknife_stats(fl_name, save_arrays=True):
    """get parameters with CI from jackknife results"""
    # set path to main folder of each jackknife result
    fdir = '{}/result/final_files/{}/'.format(root_dir, fl_name)
    # create list of filenames for saving or loading presaved data
    ftokens = ['pmf.npy', 'udl.npy', 'pi0.npy', 'clh.npy']
    # original and jackknife R^2 files
    run_id = 'fish95_euarchontoglires30_filtered'
    f_rsq = '{}/result/final_files/{}/rsq.log'.format(root_dir, run_id)
    f_jkrsq = fdir + 'jkrsq.log'

    # save data arrays generated from original files if flagged
    if save_arrays:
        pmf = []
        udl = []
        pi0 = []
        clh = []
        for fldr in os.listdir(fdir):
            f_path = fdir + fldr
            # check contents of dirs
            if os.path.isdir(f_path):
                f_list = os.listdir(f_path)
                if len(f_list):
                    c = [f for f in f_list if 'composite' in f]
                    # should only have one composite file
                    if len(c) > 1:
                        print "MULTIPLE COMPOSITES! {}".format(f_path)
                    # initialize RunStruct with composite file
                    f_init = '{}/{}'.format(f_path, c[0])
                    cst = ChromStruct('chr1', init=f_init)
                    # get dfe, udel, pi0 and CLLH from composite run
                    pmf.append(cst.uvec[0] * 1e8)
                    udl.append(sum(cst.uvec[0]) * 1e8)
                    pi0.append(cst.params[-1] / cst.fixed.tau_init)
                    clh.append(cst.stat.best_lh)

        # convert to arrays and save
        dtlists = [pmf, udl, pi0, clh]
        for ft, dt in zip(ftokens, dtlists):
            fsave = fdir + ft
            np.save(fsave, np.array(dt))

    # load data arrays
    dtfiles = [fdir + ft for ft in ftokens]
    pmf, udl, pi0, clh = [np.load(f) for f in dtfiles]

    # create the plot
    plt.figure(figsize=(12, 5))
    plt.suptitle('jackknife 2Mb windows')
    plt.subplots_adjust(left=0.08, top=0.82, wspace=1.5, right=0.9)

    # plot DFE
    plt.subplot(1, 11, (5, 8))
    plt.title('distribution of fitness effects')
    pmf_avg = np.average(pmf, axis=0)
    pmf_min = np.min(pmf, axis=0)
    pmf_max = np.max(pmf, axis=0)
    xax = np.arange(4)
    # pmf_erange = (pmf_max - pmf_min) / 2.0
    pmf_erange = np.sqrt(jackknife_variance(pmf)) / 2.0
    plt.bar(xax, pmf_avg, yerr=pmf_erange, ecolor='k')
    plt.ylabel(r'$\mu_{del} \cdot 10^8$')
    xtck = [-3.5, -3.0, -2.5, -2.0]
    plt.xticks(xax, [r'$10^{%.1f}$' % x for x in xtck])
    plt.xlabel('deleterious fitness effect')

    # plot udel
    ax_1 = plt.subplot(1, 11, 9)
    ax_1.yaxis.tick_right()
    ax_1.yaxis.set_label_position('right')
    plt.xlabel(r'$\mu_{del} \cdot 10^8$')
    # udl_erange = (udl.max()-udl.min()) / 2.0
    udl_erange = np.sqrt(jackknife_variance(udl)) / 2.0
    plt.bar(0, udl.mean(), yerr=udl_erange, ecolor='k')
    plt.xticks([])
    # plt.ylim(udl.min()-2*udl_erange, udl.max()+2*udl_erange)

    # plot pi0
    ax_2 = plt.subplot(1, 11, 10)
    ax_2.yaxis.tick_right()
    ax_2.yaxis.set_label_position('right')
    plt.xlabel(r'$\frac{\pi}{\pi_0}$')
    # pi0_erange = (pi0.max()-pi0.min()) / 2.0
    pi0_erange = np.sqrt(jackknife_variance(pi0)) / 2.0
    plt.bar(0, pi0.mean(), yerr=pi0_erange, ecolor='k')
    plt.xticks([])
    # plt.ylim(pi0.min()-2*pi0_erange, pi0.max()+2*pi0_erange)

    # plot cllh
    ax_3 = plt.subplot(1, 11, 11)
    ax_3.yaxis.tick_right()
    ax_3.yaxis.set_label_position('right')
    plt.xlabel('-composite LLH')
    clh = 1e-5 * clh
    # clh_erange = (clh.max()-clh.min()) / 2.0
    clh_erange = np.sqrt(jackknife_variance(clh)) / 2.0
    plt.bar(0, clh.mean(), yerr=clh_erange, ecolor='k')
    plt.xticks([])
    # plt.ylim(clh.min()-2*clh_erange, clh.max()+2*clh_erange)

    # plot R^2
    win, rsq = np.loadtxt(f_rsq)[:16].T
    ws = np.log10(win)
    jkrsq = np.loadtxt(f_jkrsq)[:16, 1]
    plt.subplot(1, 11, (1, 4))
    # ax_4.yaxis.tick_right()
    # ax_4.yaxis.set_label_position('right')
    # plt.plot(ws-0.05, rsq, marker='o', lw=0, label='true', ms=5)
    # plt.plot(ws+0.05, jkrsq, marker='s', lw=0, fillstyle='none', label='jackknife',
    #          markeredgewidth=1)
    plt.plot(ws-0.02, rsq, marker='s', lw=0, label='true', ms=7, alpha=0.75)
    plt.plot(ws+0.02, jkrsq, marker='s', lw=0, label='jackknife', ms=7,
             alpha=0.75)
    plt.title('variance explained')
    plt.ylabel(r'$R^2$')
    plt.xlabel('window size')
    # plot subset of x ticks
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4,5,6]])
    plt.legend(prop=dict(size=9), loc='upper left')

    # save high res figure
    fsave = fdir + 'stats_plot.png'
    plt.savefig(fsave, dpi=256)
    plt.close()
    # plt.show()


def compare_cons_subrates(filt=False):
    """compare ratio of cons/neut subrates for different % cons and n-spec"""
    # set results file directory
    pdir = root_dir + '/data/phast/hco_neutral/'

    # first count neutral subs
    nt_fmt = pdir + '{}.hco.exptotsub'
    nt_mm = mean_matrix(nt_fmt)
    # nt_bases = sum(m.sum() for m in nt_mm.values())
    # nt_sub = nt_bases - sum(m.diagonal().sum() for m in nt_mm.values())
    nt_bases = nt_mm['hg19'].sum()
    nt_sub = nt_bases - nt_mm['hg19'].diagonal().sum()
    nt_rate = 1.0 * nt_sub / nt_bases

    # create plot
    plt.figure(figsize=(12, 5))
    # plt.suptitle('conserved/neutral across conserved annotations')
    plt.subplots_adjust(left=0.08, top=0.82, right=0.9)
    xvals = np.arange(3)

    # species and percentage conserved used
    species = 'ape euarchontoglires fish'.split()
    pct_con = [0.92, 0.95, 0.98]

    # for each percentage conserved, group the ratios for different species
    for (n, pct) in enumerate(pct_con, start=1):
        p_title = '{:.0f}% conserved'.format(100*(1-pct))
        if filt:
            p_title += ' (neutral filtered)'
        pct_group = []
        for spc in species:
            # format the conserved xts filename
            cn_fmt = pdir + '{{}}.hco.{}.cons.{}.'.format(spc, pct)
            if filt:
                cn_fmt += 'filt.exptotsub'
            else:
                cn_fmt += 'exptotsub'
            # count conserved substitutions
            cn_mm = mean_matrix(cn_fmt)
            # cn_bases = sum(m.sum() for m in cn_mm.values())
            # cn_sub = cn_bases - sum(m.diagonal().sum() for m in cn_mm.values())
            cn_bases = cn_mm['hg19'].sum()
            cn_sub = cn_bases - cn_mm['hg19'].diagonal().sum()
            cn_rate = 1.0 * cn_sub / cn_bases
            pct_group.append(cn_rate / nt_rate)
        plt.subplot(1, 3, n)
        plt.title(p_title)
        plt.bar(xvals, pct_group)
        plt.xticks(xvals, [s[:4] for s in species])
        if n == 1:
            plt.ylabel('conserved/neutral substitution rates in HCO')
    sdir = root_dir + '/result/final_files/'
    fsave = sdir + 'hg19branch.cons2neut.subrate.png'
    if filt:
        fsave = sdir + 'hg19branch.cons2neut.subrate.filtered.png'
    plt.savefig(fsave, dpi=256)
    plt.close()
    # plt.show()


def annotation_overlap_venn(pct_cons, filt=False):
    pdir = root_dir + '/data/phast/hco_venn'
    cdict = defaultdict(float)
    f_fmt = '{}/{}.cons.overlap.{}.filter_{}.log'
    for ch in human_autosomes:
        f_name = f_fmt.format(pdir, ch, pct_cons, filt)
        with open(f_name, 'r') as f:
            for line in f:
                name, count = line.split()
                cdict[name] += (float(count) / 1e6)
    tot = int((cdict['ape']+cdict['euarchontoglires']+cdict['fish'])/3.0)
    lbl = 'ape euarchontoglires fish'.split()
    ape = tot - cdict['A+E'] - cdict['A+F'] + cdict['A+E+F']
    ae = cdict['A+E'] - cdict['A+E+F']
    af = cdict['A+F'] - cdict['A+E+F']
    euar = tot - cdict['A+E'] - cdict['E+F'] + cdict['A+E+F']
    ef = cdict['E+F'] - cdict['A+E+F']
    fish = tot - cdict['A+F'] - cdict['E+F'] + cdict['A+E+F']
    # sbst = (ape, euar, ae, fish, af, ef, cdict['A+E+F'])
    sbst = (ape/tot, euar/tot, ae/tot, fish/tot, af/tot, ef/tot,
            cdict['A+E+F']/tot)
    sbst = tuple(round(s, 3) for s in sbst)

    # subset order = Abc, aBc, ABc, abC, AbC, aBC, ABC
    pct_lbl = int(100-(100*pct_cons))
    titl = 'fraction overlap in top {}% conserved'.format(pct_lbl)
    # titl = 'top {}% conserved (units = Mb)'.format(int(100-(100*pct_cons)))

    # titl = 'top 8% conserved (untis = Mb)'

    if filt:
        titl += ' -- neutral filtered'
    plt.title(titl)
    v = venn.venn3(subsets=sbst, set_labels=lbl)
    # plt.show()
    fsave = '{}/venn.{}.filter_{}.png'.format(pdir, pct_cons, filt)
    plt.savefig(fsave, dpi=256)
    plt.close()


def select_utot():
    fdir = root_dir + '/result/final_files/refiltered'
    spc = 'ape euarchontoglires fish'.split()
    p_fmt = '{}/{}{}_euarchontoglires30_filtered/'
    plt.figure(figsize=(12, 5))
    plt.subplots_adjust(left=0.08, top=0.82, right=0.9)
    for (n, pct) in enumerate([0.92, 0.95, 0.98], start=1):
        pstr = int(100*pct)
        plt.subplot(1, 3, n)
        plt.title('{:.0f}% conserved'.format(100*(1-pct)))
        ymax = 0
        for (x, sp) in enumerate(spc):
            p_name = p_fmt.format(fdir, sp, pstr)
            f_list = [f for f in os.listdir(p_name) if 'composite.txt' in f]
            assert len(f_list) == 1
            f_name = p_name + f_list[0]
            utot = ChromStruct('chr1', init=f_name).stat.utot[0] * 1e8
            ymax = max(ymax, utot)
            plt.bar(x, utot, label=sp)
        plt.ylim(0, 1.25*ymax)
        if n == 1:
            plt.ylabel(r'$\mu_{del} \cdot 10^8$')
        plt.xticks(range(3), spc)
    f_save = '{}/aef_utot.png'.format(fdir)
    plt.savefig(f_save, dpi=256)
    plt.close()
    # plt.show()


def cons_dist_corr(bin_size):
    # collect files for bin size
    pth = root_dir + '/data/phast/gmap_cons_dists/'
    lbl = 'dist.bin_{:.2e}.txt'.format(bin_size)
    files = [pth + f for f in os.listdir(pth) if f.endswith(lbl)]
    # assert len(files) == 22

    # join files into single array
    arr = np.concatenate([np.loadtxt(f) for f in files]).T
    ae = pearsonr(arr[1], arr[2])
    af = pearsonr(arr[1], arr[3])
    ef = pearsonr(arr[2], arr[3])

    print 'conserved segment correlation in {:.2e} cM bins'.format(bin_size)
    print 'AE: corr={}, pval={}'.format(*ae)
    print 'AF: corr={}, pval={}'.format(*af)
    print 'EF: corr={}, pval={}\n'.format(*ef)
    return ae[0], af[0], ef[0]


def plot_fish_zeros_dist():
    f_fmt = root_dir + '/data/phast/{}.fish.zeros.txt'
    sp = 'ape primate prosimian euarchontoglires laurasiatheria mammal'.split()
    plt.figure(figsize=(6, 6))
    plt.title('distribution of phastcons scores at sites where fish = 0')
    for s in sp:
        score, count = np.loadtxt(f_fmt.format(s)).T
        count = np.cumsum(count)
        count /= count[-1]
        plt.step(score, count, label=s)
    plt.xlabel('phastcons score')
    plt.xticks(np.arange(0, 1.01, 0.1))
    plt.ylabel('cumulative fraction')
    plt.yticks(np.arange(0, 1.01, 0.1))
    plt.legend()
    f_save = root_dir + '/data/phast/fish.zeros.dist.png'
    plt.savefig(f_save, dpi=256)
    plt.close()


def figure_summary(rdir, titl):
    fdir = root_dir + '/result/final_files/{}/'.format(rdir)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    clist = [f for f in os.listdir(fdir) if ('chr1' in f) and
             (f.endswith('.txt'))]
    slist = [f for f in os.listdir(fdir) if ('predsort' in f) and
             (f.endswith('.txt'))]
    r_file = fdir + 'rsq.log'
    b_file = fdir + 'b_percentiles.txt'

    f_init = fdir + flist[0]
    f_chrom = fdir + clist[0]
    f_sort = fdir + slist[0]
    f_save = fdir[:-1] + '_summary_plot.png'
    rst = RunStruct(init=f_init)

    plt.figure(figsize=(10, 7.5))
    plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, hspace=0.3,
                        top=0.90, wspace=0.35)
    r = np.loadtxt(r_file)[13,1]
    cl = rst.stat.best_lh
    ud = rst.stat.utot[0] * 1e8
    lbl = r'$R^2=%.3f;\ -CLLH=%.3f;\ u_{del}=%.3f\cdot10^{-8}$' %(r, cl, ud)
    tit = titl + '\n'  + lbl
    plt.suptitle(tit, fontsize=16)

    # SORTED MAP 1
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    div, pi, pred = np.loadtxt(f_sort).T
    xi = np.arange(div.size)
    if 'cthresh' in rdir:
        pi0 = pred[-1]
        f_save = fdir[:-1] + '_summary_plot_v2.png'
    pi /= pi0
    pred /= pi0
    plt.subplot(3, 4, (1, 6))
    plt.axhline(y=1, color='k', ls='--', alpha=0.8)
    plt.text(50, 1.03, 'w/o background selection', ha='center', va='center')
    meanpi = pi.mean()
    plt.axhline(y=meanpi, color='darkslategray', ls=':', alpha=0.8)
    # plt.text(60, meanpi-0.03, 'mean diversity', ha='left', va='center',
    #          color='darkslategray')
    plt.text(65, meanpi - 0.03, 'mean '+r'$\pi/\pi_0$', ha='left',
             va='center', color='darkslategray')
    plt.plot(xi, pi, label='observed', color='darkslategray')
    plt.plot(xi, pred, label='predicted', color='fuchsia')
    plt.ylabel(r'$\pi/\pi_0$', fontsize=14)
    plt.ylim(0.3, 1.2)
    plt.xlabel('background selection bin')
    plt.legend(loc='lower right')

    # SORTED MAP 2
    bval = np.loadtxt(b_file)
    plt.subplot(3, 4, (3, 8))
    plt.axhline(y=1, color='k', ls='--', alpha=0.8)
    plt.text(50, 1.03, 'w/o background selection', ha='center', va='center')
    meanpi = pi.mean()
    plt.axhline(y=meanpi, color='darkslategray', ls=':', alpha=0.8)
    # plt.text(0.6, meanpi - 0.03, r'$\bar{\pi/\pi_0}$', ha='left',
    #          va='center', color='darkslategray', fontsize=16)
    plt.plot(pred, pi, label='observed', color='darkslategray')
    plt.plot(pred, pred, label='predicted', color='fuchsia')
    plt.plot([0, 1], [0, 1], label=r'$y=x$', ls='--', lw=1, color='b')
    plt.ylabel(r'$\pi/\pi_0$', fontsize=14)
    plt.ylim(0.3, 1.2)
    plt.xlabel('B value')
    plt.xlim(0.3, 1.05)
    plt.legend(loc='lower right')

    # CHROMOSOME 1 MAP
    prd, obs, num = np.loadtxt(f_chrom).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    prd /= rst.stat.meanpi
    obs /= rst.stat.meanpi
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan
    plt.subplot(3, 4, (9, 12))
    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.5)
    plt.plot(xi, prd, label='predicted', color='fuchsia', lw=1.5, alpha=0.8)
    plt.ylabel(r'$\pi/\bar{\pi}$')
    plt.ylim(0.1,1.75)
    plt.xlabel('chr1 position (Mb)')
    plt.xticks(xrange(25, 250, 50))
    plt.xlim(0, xi[-1])
    plt.legend()

    plt.savefig(f_save, dpi=256)
    plt.close()


# rdir = 'nothresh'
# titl = 'No thresholds'
# figure_summary(rdir, titl)
#
# mbsuf = '0001 001 01 05 1 15 20 25'.split()
# bnum = [b**(1/7.4) for b in 0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25]
#
# for bt in range(20, 61, 10) + range(65, 86, 5):
#     rdir = 'ape95_bth{}'.format(bt)
#     titl = r'$Optimization\ threshold:\ B \geq {}$'.format(bt / 100.0)
#     figure_summary(rdir, titl)
# for (bt, bn) in zip(mbsuf, bnum):
#     rdir = 'ape95_minbs_{}'.format(bt)
#     titl = r'$Precalculation\ threshold:\ B \geq {:.3f}$'.format(bn)
#     figure_summary(rdir, titl)
# for ct in [1] + range(5, 101, 5):
#     rdir = 'cthresh_{:02}'.format(ct)
#     titl = r'$Error\ rate\ threshold:\ c = {}$'.format(ct / 100.0)
#     figure_summary(rdir, titl)


def main():
    # flist = 'datarun_000 nff2 datarun_650 nff2_650 std_run_feb2019 ' \
    #         'nff2_bsmin010'.split()
    # llist = [r'$B \geq 0$', r'$B \geq 0\ (old)$', r'$B \geq 0.65$',
    #          r'$B \geq 0.65\ (old)$', r'$precalc\ B \geq 0.01$',
    #          r'$precalc\ B \geq 0.01\ (old)$']
    # clist = 'deeppink mediumvioletred darkturquoise darkcyan goldenrod ' \
    #         'chocolate'.split()
    # mlist = 'o o s s ^ ^'.split()
    # sname = 'compare_oldnew'
    # flist = 'datarun_000 datarun_650 std_run_feb2019 cut_dfe1 cut_dfe2'.split()
    # llist = [r'$B \geq 0$', r'$B \geq 0.65$', r'$precalc\ B \geq 0.01$',
    #          r'$t \geq 10^{-3.5}$', r'$t \geq 10^{-3}$']
    # clist = 'hotpink deeppink mediumvioletred mediumseagreen darkgreen'.split()
    # mlist = 'o o o ^ ^'.split()
    # sname = 'compare_truncated_dfe'
    # PARAMETER PLOTS
    # flist = 'std_run_feb2019 nff2_bsmin010 pr94bsmin01'.split()
    # flist = 'datarun_000 nff2 pr94bsmin00'.split()
    # llist = [r'$precalc\ B \geq 0.01\ remove\ cons\ neut$',
    #          r'$precalc\ B \geq 0.01\ keep\ con\ neuts$',
    #          r'$precalc\ B \geq 0.01\ remove\ neut\ cons$']
    # llist = ['remove cons neut', 'keep cons neut', 'remove neut cons']
    # clist = 'deeppink darkturquoise goldenrod'.split()
    # llist = [r'$B \geq 0.65\ (old)$',
    #          r'$B \geq 0.65\ (new)$',
    #          r'$B \geq 0.65\ (old\ nonCpG)$',
    #          r'$B \geq 0.65\ (new\ nonCpG)$']

    # sname = 'newdiv_vs_old_no_bound'
    # flist = 'pr94bsmin00 newdiv00 newdivCpG00'.split()
    # llist = [r'$B \geq 0.00$',
    #          r'$B \geq 0.00\ (new)$',
    #          r'$B \geq 0.00\ (new\ nonCpG)$']

    wn = 2500
    sp = 4

    # for wn in 2500, 5000, 10000:
    #     for sp in 4, 8:
    #
    #         mlist = 'o s ^ D'.split()
    #         sname = '{}sp{}'.format(sp, wn)
    #         flist = '{s}sp{w} {s}sp{w}sl {s}sp{w}cpg {s}sp{w}cpgsl'.format(s=sp, w=wn).split()
    #         llist = 'non-overlap,sliding,non-overlap nonCpG,sliding nonCpG'.split(',')
    #         clist = 'hotpink deeppink mediumvioletred purple'.split()
    #
    #         # clist = 'darkturquoise darkcyan'.split()
    #         compare_runs(flist, llist, clist, mlist, sname)
    #         # plot_observed_vs_expected_cllh(flist, llist, clist, sname)
    #         # plot_obex_regression('/Users/davidmurphy/GoogleDrive/linked_selection/'
    #     #                      'result/final_files/pr94bsmin00/xparam.obex_0.65.npy')

    # mlist = 'o s ^ D > < ^ *'.split()
    # sname = 'phylo_depth'
    # flist = 'ape93 pri94 prosim95 euar95 laur95 mamm95 fish95'.split()
    # llist = '5-ape 4.5%,8-primate 4.1%,12-prosim 4.0%,25-euarchon. 4.3%,' \
    #         '50-laurasia. 4.3%,61-mammal 4.4%,98-fish 4.4%'.split(',')
    # clist = 'salmon indianred firebrick darkturquoise teal goldenrod ' \
    #         'royalblue'.split()
    npct = '30'
    cpct = '95'
    mlist = 'o s ^ D > < ^ *'.split()
    sname = 'vary_cons{}_neut{}_filtered'.format(cpct, npct)
    # sname = 'compare_fish_ape_95'
    # spc = 'ape fish'.split()
    # flist = []
    # llist = []
    spc = 'ape primate prosimian euarchontoglires laurasiatheria mammal ' \
             'fish'.split()
    # for s in spc:
    #     f1 = '{}{}_euarchontoglires{}'.format(s, cpct, npct)
    #     f2 = '{}{}_euarchontoglires{}_filtered'.format(s, cpct, npct)
    #     flist.append(f1)
    #     flist.append(f2)
    #     llist.append('{} {}%'.format(s, cpct))
    #     llist.append('f {} {}%'.format(s, cpct))

    flist = ['{}{}_euarchontoglires{}_filtered'.format(s, cpct, npct) for s in spc]
    llist = ['{} {}%'.format(s, cpct) for s in spc]
    clist = 'salmon indianred firebrick darkturquoise teal goldenrod ' \
            'royalblue'.split()
    compare_runs(flist, llist, clist, mlist, sname)
    # individual_rsq_plot(flist, llist, sname)

    # clist = 'firebrick goldenrod royalblue'.split()
    # cons_2 = [str(n) for n in range(15, 96, 10) + [100]]
    # for cn in cons_2:
    #     flist = ['fish{}_euarchontoglires{}'.format(p, cn) for p in 92, 95, 98]
    #     sname = '{} percent neutral cutoff'.format(cn)
    #     llist = ['fish {}%'.format(p) for p in 92, 95, 98]
    #     compare_runs(flist, llist, clist, mlist, sname)
    # clist = 'darkturquoise darkcyan'.split()
    # for pct in xrange(20, 101, 10):
    #     fl1 = ['{}92_euarchontoglires{}'.format(s, pct) for s in spc]
    # compare_runs(flist, llist, clist, mlist, sname)
    # individual_rsq_plot(flist, llist, sname)
    # grouped_rsq()
    # grouped_rsq_vary_neut()

    # # set intervals used for estimating substitution rates
    # intvl_1 = [(0, 30), (30, 40), (40, 50), (50, 60), (60, 70),
    #            (70, 80), (80, 90), (90, 100)]
    # intvl_2 = [(0, x) for x in xrange(20, 101, 10)]
    # plot_phylofit_neut_subrates(intvl_1, 'neut_subrates_sum_all_1')
    # plot_phylofit_neut_subrates(intvl_2, 'neut_subrates_sum_all_2')
    # for pc in [92, 95, 98]:
    #     for filt in [True, False]:
    #         grouped_rsq_vary_spc_cons_neut(pc, filt)


clist = 'darkorange purple gold fuchsia limegreen blue firebrick'.split()
mlist = 'o s ^ D > < ^ * H'.split()
# sname = 'exoncons annotations'
# flist = ['ape95_minbs_01', 'ape_exnex']
# llist = ['ape 5%', 'conserved exonic']
# clist = ['fuchsia', 'orange', 'purple']
# mlist = 'o s'.split()
spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()
flist = ['{}_cons95_clean'.format(sp) for sp in spec]
llist = spec
sname = 'clean_run'
compare_runs(flist, llist, clist, mlist, sname)
# individual_rsq_plot(flist, llist, sname)

