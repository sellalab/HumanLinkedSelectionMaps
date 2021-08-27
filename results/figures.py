from classes.runstruct import RunStruct, root_dir, goog_dir, izip, np, os, re
from results.moving_avg import moving_avg
from scipy.stats import pearsonr as pcor
# from functions import swap_root
import matplotlib.pyplot as plt
import matplotlib.colors as clr

import seaborn


__author__ = 'davidmurphy'


result_dir = root_dir + '/result'
figure_dir = goog_dir + '/linked_selection/murphy_et_al/figures'
cols = ['HotPink', 'MediumVioletRed', 'Coral', 'OrangeRed', 'Cyan', 'Teal',
        'MediumAquamarine', 'DarkGreen', 'Indigo', 'MediumOrchid']
mrks = 'h > D v o ^ s < p x'.split()


def sorted_plot():

    rdir = root_dir + '/result/final_files/'
    # folder = 'alternateclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    # folder = 'popsclean'
    folder = 'sensitivity_tests'
    d = rdir + folder

    files = ['{}/{}'.format(d, f) for f in os.listdir(d) if f.endswith('txt')]

    # titles = 'CADD McVicker Genic/Nongenic-cons'.split()
    # titles = [f.split('/')[-1].split('.')[0] for f in files]
    # titles = ['{}%'.format(p) for p in xrange(90, 100)]
    titles = ['conserved', 'composite ex/nex', 'conserved-high']

    num_bins = 500

    for (f, t) in zip(files, titles):
        plot_sorted_map(init_file=f, folder=folder, title=t, num_bins=num_bins)
    plt.show()


def plot_sorted_map(init_file, folder, title, num_bins):
    from functions import swap_root
    swap_root(init_file)
    rst = RunStruct(init=init_file)
    arrtmp = '{}/result/sort_map/{}/{}.{}.sorted.{}bins.npy'
    arr_file = arrtmp.format(rst.root, folder, rst.neut, rst.label, num_bins)
    figdir = '{}/linked_selection/murphy_et_al/figures/sort_{}'.format(goog_dir,
                                                                       folder)
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    figtmp = '{}/{}.{}.sortmap.{}bins.png'
    fig_file = figtmp.format(figdir, rst.neut, rst.label, num_bins)

    # rst.label = '{}.BS{}.{}.CS{}.{}'.format(rst.token, rst.bnum, rst.bsgrid, rst.cnum, rst.csgrid)
    # arr_file = '{}/result/sort_map/primate90to99/{}.{}.sorted.{}bins.npy'.format(rst.root, rst.neut, rst.label,
    #                                                                              num_bins)
    # arr = basic_sort(rst=rst, num=num_bins)
    # np.save(arr_file, arr)
    arr = np.load(arr_file)

    # dv, pi, pr, ar1, ar2, hnd = arr.T
    dv, pi, pr = arr.T

    meanpi = rst.stat.meanpi / rst.stat.meanpi0
    x = np.arange(len(dv))
    pi /= rst.stat.meanpi0
    pr /= rst.stat.meanpi0
    # ar1 /= ar1.mean()
    # ar1 *= meanpi
    # ar2 /= ar2.mean()
    # ar2 *= meanpi
    # hnd /= hnd.mean()
    # hnd *= meanpi

    corr, pval = pcor(pi, pr)

    gray = 'DarkSlateGray'
    figwidth = 10
    fig = plt.figure(figsize=(figwidth, 7))
    plt.subplots_adjust(left=0.13, bottom=0.12, right=0.78, top=0.90)
    fig.add_subplot(111)

    plt.plot(x, pi, label='observed diversity', color=gray, lw=3)
    plt.plot(x, pr, label='predicted', color='Fuchsia', alpha=0.75, lw=3)

    # plt.plot(x, ar1, label='archaic haplotype probability (Song)', color='DodgerBlue', alpha=0.75, lw=3)
    # plt.plot(x, ar2, label='archaic haplotype probability (Reich)', color='Purple')
    # plt.plot(x, hnd, label='scaled YRI-Altai divergence', color='FireBrick', alpha=0.75, lw=3)

    message = ' {}\n {:>13}'.format(*'no background;selection'.split(';'))
    plt.axhline(y=1, color='k', ls='--', lw=2.0)
    # plt.text((0.5 * num_bins), 1.005, 'w/o background selection', ha='center', va='bottom', fontsize=18)
    plt.text((1.025 * num_bins), 1, message, ha='left', va='center', fontsize=20)

    message = '   {:>6}\n   {}'.format(*'mean;diversity'.split(';'))
    plt.axhline(y=meanpi, color=gray, ls=':', lw=2.0)
    # plt.text((0.8 * num_bins), meanpi - 0.005, 'mean diversity', ha='center', va='top', fontsize=18, color=gray)
    plt.text((1.025 * num_bins), meanpi, message, ha='left', va='center', fontsize=20, color=gray)

    # Pearson correlation
    plt.text((0.5 * num_bins), 0.55, r'$Pearson\ R^2 = {:.4f}$'.format(corr), ha='center', va='center', fontsize=20)

    # title w/ species or cons info
    plt.title(title, fontsize=20)
    # plt.title(rst.token, fontsize=20)
    # plt.title(init_file.split('/')[-1].split('.')[0], fontsize=20)

    plt.xlabel('strength of background selection', fontsize=20)
    plt.xticks(fontsize=18)
    plt.xlim((-0.01 * num_bins), (1.01 * num_bins))

    plt.ylabel('scaled diversity', fontsize=20, labelpad=10)
    plt.yticks(fontsize=18)
    plt.ylim(0.37, 1.33)

    plt.legend(prop={'size': 18}, loc='upper left', ncol=1, numpoints=3, borderaxespad=0.1)

    # plt.savefig(fig_file, dpi=2560/figwidth)
    # plt.show()


def collated_plot(focal):
    """collated plots around cds"""

    rdir = root_dir + '/result/final_files/'
    # folder = 'pr90to99clean'
    # folder = 'phyloclean'
    # folder = 'alternateclean'
    # folder = 'popsclean'
    folder = 'sensitivity_tests/100bins'

    d = rdir + folder
    files = ['{}/{}'.format(d, f) for f in os.listdir(d) if f.endswith('txt')]
    _ = [swap_root(f) for f in files]
    rs = [RunStruct(init=f) for f in files]

    collated_dir = '{}/result/collate/{}'.format(rs[0].root, folder)
    coltmp = '{dd}/{neut}.{label}.{fc}.collated.npy'
    colfiles = [coltmp.format(dd=collated_dir, fc=focal, **r.dict) for r in rs]
    cx = [np.load(f).T for f in colfiles]

    labels = ['conserved', 'composite ex/nex', 'conserved-high']
    # labels = [f.split('/')[-1].split('.')[0] for f in files]
    # labels = ['{}%'.format(p) for p in xrange(1, 11)]
    # labels = 'CADD McVicker 08-primate'.split()

    focal_dict = dict(syn='synonymous substitution',
                      nonsyn='nonsynonymous substitution',
                      primate_conserved='primate conserved substitution',
                      cds_midpoints='protein coding segment',
                      nc_cons='conserved noncoding segment')

    # plt.figure(figsize=(10, 7))
    # plt.subplots_adjust(left=0.12, bottom=0.12, right=1.0, top=1.0, wspace=0, hspace=0)
    # bins, div, pi, pr, cnts = cx[0]
    # obs = pi / (cnts * rs[0].stat.meanpi)
    # plt.plot(bins, obs, color='gray', lw=2)

    i = 0
    for (r, c, l) in zip(rs, cx, labels):
        plt.figure(figsize=(10, 7))
        plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=1.0)
        bins, div, pi, pr, cnts = c
        obs = pi / (cnts * r.stat.meanpi0)
        pred = pr / (cnts * r.stat.meanpi0)
        ymin = min(pred.min(), obs.min())
        ymax = max(pred.max(), obs.max())

        plt.plot(bins, obs, color='gray', lw=2)
        plt.plot(bins, pred, color=cols[i], alpha=0.75, lw=1.5)
        i += 1

        plt.xticks(fontsize=20)
        plt.xlabel('cM to nearest {}'.format(focal_dict[focal]), fontsize=22,
                   labelpad=10)
        plt.xlim(-rs[0].vars.collated_span, rs[0].vars.collated_span)

        plt.yticks(fontsize=20)
        plt.ylabel('diveristy', fontsize=22, labelpad=10)
        yspan = ymax - ymin
        # plt.ylim(ymin - (0.1 * yspan), ymax + (0.2 * yspan))

        legend_text = ['Observed diversity', l]
        plt.legend(legend_text, prop={'size': 16}, ncol=3, loc='upper center')

        figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/collated_{}'.format(folder)
        if not os.path.isdir(figdir):
            os.mkdir(figdir)
        fig_file = '{}/{}.{}.png'.format(figdir, l, focal)
        # plt.savefig(fig_file, dpi=256)
        # plt.show()
    plt.show()


def rsq_plot():
    rdir = root_dir + '/result/rsq_maps/'
    # folder = 'alternateclean'
    # folder = 'partitions'
    # folder = 'variants'
    # folder = 'popsclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    folder = 'sensitivity_tests'
    d = rdir + folder
    files = ['{}/{}'.format(d, f) for f in os.listdir(d) if f.endswith('txt')]
    # files = ['{}/{}'.format(d, f) for f in os.listdir(d)][::-1]

    # labels = ['{}%'.format(p) for p in xrange(1, 11)]
    # labels = [f.split('/')[-1].split('.')[0] for f in files]
    # labels = ['CEU-cons', 'CEU-cons-coding/noncoding', 'LWK-cons',
    # 'LWK-cons-coding/noncoding', 'YRI-cons', 'YRI-cons-coding/noncoding']
    # labels = ['YRI-CADD', 'YRI-cons', 'YRI-cons-coding/noncoding']
    # labels = ['conserved', 'conserved-1K-fix', 'conserved-1K',
    #           'composite', 'composite-1K-fix', 'composite-1K']
    labels = ['conserved', 'composite ex/nex', 'conserved-high']

    # arrays = [np.loadtxt(f)[:-2].T for f in files]
    arrays = [np.loadtxt(f)[:-2].T for f in files]  # REVERSED!
    win = np.log10(arrays[0][0])

    plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.09, bottom=0.1, right=1.0, top=1.0)

    space = 0.013
    shift = -space * 0.5 * len(arrays)
    # shift = 0

    for (i, arr) in enumerate(arrays):
        facecol = clr.colorConverter.to_rgba(cols[i], alpha=0.25)
        edgecol = clr.colorConverter.to_rgba(cols[i], alpha=0.75)
        mdic = dict(markerfacecolor=facecol, markeredgecolor=edgecol,
                    markeredgewidth=1.5, color=edgecol, lw=0.0)
        # mk = mrks[i]
        mk = 'o'

        plt.plot(win + shift, arr[1], marker=mk, ms=7.5, **mdic)
        shift += space

    partitions = win[:-1] + (win[1:] - win[:-1]) / 2
    for x in partitions:
        plt.axvline(x=x, ls='--', lw='0.8', color='k', alpha=0.75)

    plt.yticks(fontsize=20)
    plt.ylabel('variance explained', fontsize=22, labelpad=10)

    xlabels = ['10kb'] + ['']*6 + ['100kb'] + ['']*5 + ['1Mb', '', '']
    plt.xticks(win, xlabels, fontsize=20)
    plt.xlabel('window size', fontsize=22, x=0.474)
    plt.xlim(win[0] - (win[1] - win[0]) / 2, win[-1] + (win[-1] - win[-2]) / 2)

    leg = plt.legend(labels, prop={'size': 16}, frameon=True, ncol=1)
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_alpha(1.0)

    figdir = goog_dir + '/linked_selection/murphy_et_al/figures'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/variance_explained/{}.png'.format(figdir, folder)
    # plt.savefig(fig_file, dpi=213)
    plt.show()


def dfe_plot():
    rdir = root_dir + '/result/final_files/'
    # folder = 'alternateclean'
    # folder = 'popsclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    # folder = 'genicpartclean'
    folder = 'sensitivity_tests'
    d = rdir + folder
    files = ['{}/{}'.format(d, f) for f in os.listdir(d) if f.endswith('txt')]
    # files = ['{}/{}'.format(d, f) for f in os.listdir(d)][::-1]

    rs = [RunStruct(init=f) for f in files][::-1]

    labels = ['conserved', 'composite ex/nex', 'conserved-high']
    # labels = ['conserved', 'conserved-1K-fix', 'conserved-1K',
    #           'composite', 'composite-1K-fix', 'composite-1K']
    # labels = ['{}%'.format(p) for p in xrange(1, 11)]
    # labels = [f.split('/')[-1].split('.')[0] for f in files]
    # labels = [['CEU genic cons', 'LWK genic cons', 'YRI genic cons'],
    #           ['CEU nongenic cons', 'LWK nongenic cons', 'YRI nongenic cons']]

    figwidth = 10
    plt.figure(figsize=(figwidth, 7))
    plt.subplots_adjust(left=0.1, bottom=0.13, right=1, top=1, wspace=0.1,
                        hspace=0.12)

    x = np.arange(rs[0].bsgrid)
    for i in xrange(rs[0].bnum):
        plt.subplot(rs[0].bnum, 1, i+1)
        arr = np.array([r.stat.upmf[i] for r in rs])
        # x = np.arange(len(arr[0]))
        width = 0.8 / len(arr)
        shift = -(0.5 * len(arr) * width)
        for (j, row) in enumerate(arr):
            # plt.plot(x, np.cumsum(row), '-s', color=cols[i])
            plt.bar(x + shift, row, width=width, align='edge', color=cols[j])
            shift += width

        plt.xticks(x, ['']*len(x))
        plt.yticks(fontsize=20)

        lab = labels[i] if rs[0].bnum > 1 else labels
        leg = plt.legend(lab, prop={'size': 16}, frameon=True, loc='upper left')
        frame = leg.get_frame()
        frame.set_facecolor('white')
        frame.set_alpha(1.0)

    yshift = 1 if rs[0].bnum > 1 else 0.5
    plt.ylabel('proportion of deleterious mutations', fontsize=22, y=yshift)
    if folder == 'phyloclean':
        plt.ylim(0, 0.65)
    xstr = r'$\mathrm{10^{%.1f}}$'
    plt.xticks(x, [xstr % s for s in np.arange(-4.5, -1.5, 0.5)], fontsize=20)
    plt.xlabel('deleterious mutation fitness effects', fontsize=22, labelpad=15)

    leg = plt.legend(labels, prop={'size': 16}, frameon=True, loc='upper left')
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_alpha(1.0)

    figdir = goog_dir + '/linked_selection/murphy_et_al/figures/dfe'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, folder)
    # plt.savefig(fig_file, dpi=2560/figwidth)
    plt.show()


def udel_plot():
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/final_files/'
    # folder = 'alternateclean'
    # folder = 'popsclean'
    folder = 'phyloclean'
    # folder = 'pr90to99clean'
    # folder = 'genicpartclean'

    d = rdir + folder
    files = ['{}/{}'.format(d, f) for f in os.listdir(d)]
    # files = ['{}/{}'.format(d, f) for f in os.listdir(d)][::-1]
    rs = [RunStruct(init=f) for f in files]
    print '\n'.join(r.label for r in rs)

    # labels = [['{}%'.format(p) for p in xrange(1, 11)]]
    # labels = [[r.neut for r in rs]]
    # labels = [[f.split('/')[-1].split('.')[0] for f in files]]
    labels = [['CEU genic cons', 'LWK genic cons', 'YRI genic cons'],
              ['CEU nongenic cons', 'LWK nongenic cons', 'YRI nongenic cons']]

    bnum = rs[0].bnum

    # DELETERIOUS MUTATION RATES
    plt.figure(figsize=(10, 7))
    for j in xrange(bnum):
        plt.subplot(bnum, 1, j+1)
        plt.subplots_adjust(left=0.12, bottom=0.04, right=1, top=0.95, wspace=0.1, hspace=0.12)
        for (i, r) in enumerate(rs):
            plt.bar(left=i, height=r.stat.utot[j]*1e8, color=cols[i])

        if folder == 'popsclean':
            ncol = 4
            loc = 'upper center'
        elif folder == 'pr90to99clean':
            ncol = 2
            loc = 'upper right'
        elif folder == 'genicpartclean':
            ncol = 1
            loc = 'best'
            plt.ylim(0, [1.2, 3][j])
        else:
            ncol = 3
            loc = 'upper center'
            plt.ylim(0, 2.3)

        plt.xticks([])
        plt.yticks(fontsize=20)
        leg = plt.legend(labels[j], prop={'size': 16}, frameon=True, loc=loc, ncol=ncol)
        frame = leg.get_frame()
        frame.set_facecolor('white')
        frame.set_alpha(1.0)
    plt.ylabel(r'$\mathrm{deleterious\ \mu\ (x10^{-8})}$', fontsize=22, y=0.5*bnum)
    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/udel'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, folder)
    plt.savefig(fig_file, dpi=256)

    # MEAN FITNESS EFFECT OF DELETERIOUS MUTATIONS
    plt.figure(figsize=(10, 7))
    pad = 0
    for j in xrange(bnum):
        plt.subplot(bnum, 1, j+1)
        plt.subplots_adjust(left=0.12, bottom=0.04, right=1, top=0.95, wspace=0.1, hspace=0.12)
        for (i, r) in enumerate(rs):
            plt.bar(left=i, height=r.stat.tmean[j]*1e3, color=cols[i])

        if folder == 'popsclean':
            ncol = 2
            loc = 'upper right'
        elif folder == 'pr90to99clean':
            ncol = 2
            loc = 'upper right'
        elif folder == 'genicpartclean':
            ncol = 1
            loc = 'upper left'
            plt.ylim(0, [0.45, 3.5][j])
            pad = 10
        else:
            ncol = 3
            loc = 'upper center'
            plt.ylim(0, 7)

        plt.xticks([])
        plt.yticks(fontsize=20)
        leg = plt.legend(labels[j], prop={'size': 16}, frameon=True, loc=loc, ncol=ncol)
        frame = leg.get_frame()
        frame.set_facecolor('white')
        frame.set_alpha(1.0)
    plt.ylabel(r'$\mathrm{mean\ fitness\ effect\ (x10^{-3})}$', fontsize=22, y=0.5*bnum, labelpad=pad)
    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/mean_fitnes_effects'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, folder)
    plt.savefig(fig_file, dpi=256)

    # OVERALL REDUCTION IN DIVERSITY
    lab = labels
    plt.figure(figsize=(10, 7))
    plt.subplots_adjust(left=0.12, bottom=0.04, right=1, top=0.95, wspace=0.1, hspace=0.12)
    for (i, r) in enumerate(rs):
        plt.bar(left=i, height=r.stat.diversity_reduction, color=cols[i])

    if folder == 'popsclean':
        ncol = 4
        loc = 'upper center'
    elif folder == 'pr90to99clean':
        ncol = 5
        loc = 'upper center'
    elif folder == 'genicpartclean':
        ncol = 1
        loc = 'upper left'
        lab = [r.neut for r in rs]
        plt.ylim(0, 0.22)
    else:
        ncol = 3
        loc = 'upper center'

    plt.ylabel('diversity reduction', fontsize=22)
    plt.xticks([])
    plt.yticks(fontsize=20)
    leg = plt.legend(lab, prop={'size': 16}, frameon=True, loc=loc, ncol=ncol)
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_alpha(1.0)
    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/diversity_reduction'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, folder)
    # plt.savefig(fig_file, dpi=256)
    plt.show()


def cons_filters(sp):
    dr = '{}/cons_filters/{}/'.format(result_dir, sp)
    fs = [dr + f for f in os.listdir(dr) if f.startswith(sp)]
    ar = [np.loadtxt(f, dtype=str).T for f in fs]
    rp = [re.search('\.(\d\d)\.txt$', f).group(1) for f in fs]
    nm = ['{}% threshold'.format(p) for p in rp]

    width = 0.9 / len(ar)
    shift = -0.4
    x = np.arange(len(ar[0][0])-1)
    ymin = 1.0
    ymax = 1.0

    plt.figure(figsize=(13, 4))
    plt.subplots_adjust(left=0.06, bottom=0.08, right=1, top=0.97, wspace=0.1,
                        hspace=0.12)
    for (p, a) in zip(nm, ar):
        a = a[1].astype(float)
        frac = a[1:]/a[0]
        ymin = min(ymin, frac.min())
        ymax = max(ymax, frac.max())
        plt.bar(left=x+shift, height=frac, width=width, align='edge', label=p)
        shift += width

    plt.ylabel('fraction remaining', fontsize=12)
    plt.yticks(fontsize=12)
    yspan = ymax - ymin
    plt.ylim(ymin - 0.1*yspan, 1.0)

    # plt.xlabel('phylogeny', fontsize=22)
    plt.xticks(x, ar[0][0][1:], fontsize=12)

    plt.legend(ncol=4, loc='upper center', prop={'size': 12})
    fimg = figure_dir + '/filter_intersections/ape80to95pct.png'
    plt.savefig(fimg, dpi=256)
    # plt.show()


def alignment_coverage():
    fname = result_dir + '/100way_species_coverage/autosomes.100way.cov.txt'
    with open(fname, 'r') as f:
        a = [l.split() for l in f.readlines()]
        s, n = [], []
        for l in a:
            s.append(' '.join(l[:-1]))
            n.append(float(l[-1]))
        n = np.array(n)
    auto_len = 2881033286.0
    plt.figure(figsize=(15, 10))
    plt.subplots_adjust(right=1, left=0.04, bottom=0.18, top=0.97)
    plt.bar(left=range(100), height=n / auto_len, color='seagreen')
    plt.xticks(range(100), s, rotation=90, fontsize=10)
    plt.yticks(np.arange(0, 1.01, 0.05), rotation=45, fontsize=10)
    plt.ylabel('autosome coverage', fontsize=10)
    plt.xlim(-1, 100)
    plt.ylim(0, 0.95)
    # plt.axhline(y=0.5, ls='--', lw=1.0)
    fout = '/Users/davidmurphy/Desktop/autosomes.cov.png'
    plt.savefig(fout, dpi=256)
    # plt.show()


def bvalue_distributions(nbins, anno):
    dname = '/result/b_distributions/{}bin_low_udel'.format(nbins)
    d = root_dir + dname
    files = ['{}/{}'.format(d, f) for f in os.listdir(d) if anno in f]

    plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.08, bottom=0.11, right=0.96, top=0.95)

    for (i, f) in enumerate(files):
        t = re.search('t(0\.0\d{7})\.txt$', f).group(1)
        b, n = np.loadtxt(f).T
        n /= n.sum()
        l = r'$\mathrm{10^{%.1f}}$' % np.log10(float(t))
        plt.hist(b, bins=np.arange(0, nbins+1), weights=n, label=l,
                 histtype='step', lw=1, cumulative=1, alpha=True,
                 color=cols[i])

    tit = 'distribution of B values over {} bins in {}'.format(nbins, anno)
    plt.title(tit, fontsize=18)
    plt.xlabel('B values', fontsize=18)
    plt.xticks(fontsize=16)
    plt.ylabel('proportion of B values', fontsize=18)
    plt.yticks(fontsize=16)
    leg = plt.legend(loc='upper left', prop=dict(size=16),
                     title='selection coefficients')
    leg.get_title().set_fontsize('16')

    plt.show()


def collated_plot_2(fcols, finits, focal, cols, labs, fimg=None):
    """collated plot with paths and labels specificed in args"""
    for f in finits:
        swap_root(f)
    rs = [RunStruct(init=f) for f in finits]
    cs = [np.load(f) for f in fcols]
    # OLD ARRAY FORMATTING:
    # bins, pi, div, pr, cnts = np.load(fcol).T

    plt.figure(figsize=(10, 7))
    plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=0.95)
    # plt.title('divide by mean pi from pr95 run')
    # plt.title('old data')

    nopi = True
    meanpi = rs[0].stat.meanpi
    meanpi0 = rs[0].stat.meanpi0
    for (rst, cx, cl, lb) in izip(rs, cs, cols, labs):

        # bins, div, pi, pr, cnts = cx.T
        bins, pi, div, pr, cnts = cx.T

        # p0 = 0.00134677376667
        # obs = pi / np.multiply(cnts, p0)
        # pred = pr / np.multiply(cnts, p0)

        # obs = pi / (cnts * rst.stat.meanpi)
        # pred = pr / (cnts * rst.stat.meanpi)

        obs = pi / (cnts * meanpi0)
        pred = pr / (cnts * meanpi0)

        if nopi:
            label = 'Observed'
            plt.plot(bins, obs, color='darkslategray', lw=2, label=label)
            nopi = False

        plt.plot(bins, pred, color=cl, alpha=0.75, lw=2, label=lb)



        # ymin = min(pred.min(), obs.min())
        # ymax = max(pred.max(), obs.max())
        # yspan = ymax - ymin
        # plt.ylim(ymin - (0.1 * yspan), ymax + (0.2 * yspan))

    plt.xticks(fontsize=22)
    xlabel = 'distance to nearest {} (cM)'.format(focal)
    plt.xlabel(xlabel, fontsize=24, labelpad=10)
    plt.xlim(-rst.vars.collated_span, rst.vars.collated_span)

    plt.yticks(fontsize=22)
    plt.ylabel('scaled diversity', fontsize=24, labelpad=10)
    plt.ylim(0.61, 1)
    plt.legend(prop={'size': 22}, ncol=1, loc='lower left')

    # optional save
    if fimg:
        plt.savefig(fimg, dpi=256)

    plt.show()


def special_rsq_plot():
    xlabels = '5%-ape (4),5%-primate (8),5%-prosimian (12),' \
              '5%-supraprimate (25),5%-laurasiatheria (50),' \
              '5%-afrotheria (57),5%-mammal (61),5%-birds (75)'.split(',') + \
              ['{}%-primate (8)'.format(x) for x in xrange(1, 11)]
    rsq_values = [0.497945, 0.496878, 0.503170, 0.509512, 0.507454, 0.507095,
                  0.507072, 0.506247, 0.485339, 0.490420, 0.494193, 0.495884,
                  0.496878, 0.494752, 0.493670, 0.491966, 0.489348, 0.486947]
    x = np.arange(len(rsq_values))
    plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.12, bottom=0.47, right=0.95, top=0.95,
                        wspace=0.06, )
    col1, col2 = 'darkslategray', 'darkturquoise'

    facecol = clr.colorConverter.to_rgba(col1, alpha=0.25)
    edgecol = clr.colorConverter.to_rgba(col1, alpha=0.75)
    mdic = dict(markerfacecolor=facecol, markeredgecolor=edgecol,
                markeredgewidth=1.5, color=edgecol, lw=0.0)
    plt.subplot(121)
    plt.plot(x[:8], rsq_values[:8], marker='o', ms=15, **mdic)
    plt.xticks(x[:8], xlabels[:8], fontsize=22, rotation=90)
    plt.yticks(fontsize=22)
    plt.ylabel('variance explained', fontsize=24, labelpad=10)
    plt.ylim(0.482, 0.512)

    plt.subplot(122)
    facecol = clr.colorConverter.to_rgba(col2, alpha=0.25)
    edgecol = clr.colorConverter.to_rgba(col2, alpha=0.75)
    mdic = dict(markerfacecolor=facecol, markeredgecolor=edgecol,
                markeredgewidth=1.5, color=edgecol, lw=0.0)
    plt.plot(x[8:], rsq_values[8:], marker='o', ms=15, **mdic)
    plt.xticks(x[8:], xlabels[8:], fontsize=22, rotation=90)
    plt.yticks(color='white', alpha=0)
    plt.ylim(0.482, 0.512)
    plt.show()


def main():
    # sorted_plot()
    rsq_plot()
    # dfe_plot()
    # udel_plot()
    # for s in 'syn nonsyn cds_midpoints'.split():
    #     collated_plot(focal=s)
    # collated_plot('nonsyn')
    # cons_filters('ape')
    # alignment_coverage()
    # bvalue_distributions(100, anno='cons95_Seg')

    # cp = '{}/result/collate/'.format(root_dir)
    # fcols = [cp + 'YRI.pr95.cleanrun.BS1.6.CS0.0.nonsyn.collated.npy',
    #          cp + 'YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy']
    cp = '{}/result/collate/old_collated_data/nsSub/'.format(root_dir)
    # fcols = [cp + 'filter.only.keep.BGC.BS1CS0.nsSub.YRI.npy',
    #          cp + 'mcvicker.map.BS1CS0.nsSub.YRI.npy']
    # ipath = '{}/result/final_files/alternateclean/'.format(root_dir)
    # finits = [ipath + 'YRI.pr95.cleanrun.BS1.6.CS0.0.170901163420.final.txt',
    #           ipath + 'YRI.mcvicker.ref.clean.BS1.1.CS0.0.reformat.final.txt']
    # colors = ['fuchsia', 'darkturquoise']
    # labels = ['Our prediction', 'McVicker\'s prediction']
    # collated_plot_2(fcols, finits, 'amino acid substitution', colors, labels)

    fcols = [cp + 'mcvicker.map.BS1CS0.nsSub.YRI.npy']
    ipath = '{}/result/final_files/alternateclean/'.format(root_dir)
    finits = [ipath + 'YRI.mcvicker.ref.clean.BS1.1.CS0.0.reformat.final.txt']
    colors = ['darkturquoise']
    labels = ['McVicker\'s prediction']
    # collated_plot_2(fcols, finits, 'amino acid substitution', colors, labels)
    # special_rsq_plot()

if __name__ == '__main__':
    main()
