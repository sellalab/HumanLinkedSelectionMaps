from results.moving_avg import moving_avg
from classes.mapstruct import RunStruct, MapStruct, np, os
from scipy.stats import pearsonr as pcor
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import seaborn

__author__ = 'davidmurphy'

result_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/'
cols = 'HotPink MediumVioletRed Coral OrangeRed Cyan Teal MediumAquamarine DarkGreen Indigo MediumOrchid'.split()
mrks = 'h > D v o ^ s < p x'.split()


def sorted_plot():

    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/final_files/'
    folder = 'alternateclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    # folder = 'popsclean'
    d = rdir + folder

    files = ['{}/{}'.format(d, f) for f in os.listdir(d)]

    titles = 'CADD McVicker Genic/Nongenic-cons'.split()
    # titles = [f.split('/')[-1].split('.')[0] for f in files]
    # titles = ['{}%'.format(p) for p in xrange(90, 100)]

    num_bins = 500

    for (f, t) in zip(files, titles):
        plot_sorted_map(init_file=f, folder=folder, title=t, num_bins=num_bins)
    # plt.show()


def plot_sorted_map(init_file, folder, title, num_bins):
    from functions import swap_root
    swap_root(init_file)
    rst = RunStruct(init=init_file)
    arr_file = '{}/result/sort_map/{}/{}.{}.sorted.{}bins.npy'.format(rst.root, folder, rst.neut, rst.label, num_bins)
    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/sort_{}'.format(folder)
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.{}.sortmap.{}bins.png'.format(figdir, rst.neut, rst.label, num_bins)

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
    plt.subplots_adjust(left=0.13, bottom=0.12, right=0.78, top=0.90, wspace=0.05, hspace=0.2)
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

    plt.savefig(fig_file, dpi=2560/figwidth)


def collated_plot(focal):
    """collated plots around cds"""

    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/final_files/'
    # folder = 'pr90to99clean'
    # folder = 'phyloclean'
    # folder = 'alternateclean'
    folder = 'popsclean'

    init_dir = rdir + folder
    init_files = ['{}/{}'.format(init_dir, f) for f in os.listdir(init_dir)]
    # init_files = ['{}/{}'.format(init_dir, f) for f in os.listdir(init_dir)][::-1]
    rs = [RunStruct(init=f) for f in init_files]

    collated_dir = '{}/result/collate/{}'.format(rs[0].root, folder)
    collated_files = ['{dd}/{neut}.{label}.{fc}.collated.npy'.format(dd=collated_dir, fc=focal, **r.dict) for r in rs]
    cx = [np.load(f).T for f in collated_files]

    labels = [f.split('/')[-1].split('.')[0] for f in init_files]
    # labels = ['{}%'.format(p) for p in xrange(1, 11)]
    # labels = 'CADD McVicker 08-primate'.split()

    focal_dict = dict(syn='synonymous substitution', nonsyn='nonsynonymous substitution',
                      primate_conserved='primate conserved substitution', cds_midpoints='protein coding segment',
                      nc_cons='conserved noncoding segment')

    # plt.figure(figsize=(10, 7))
    # plt.subplots_adjust(left=0.12, bottom=0.12, right=1.0, top=1.0, wspace=0, hspace=0)
    # bins, div, pi, pr, cnts = cx[0]
    # obs = pi / (cnts * rs[0].stat.meanpi)
    # plt.plot(bins, obs, color='gray', lw=2)

    i = 0
    for (r, c, l) in zip(rs, cx, labels):
        plt.figure(figsize=(10, 7))
        plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=1.0, wspace=0, hspace=0)
        bins, div, pi, pr, cnts = c
        obs = pi / (cnts * r.stat.meanpi0)
        pred = pr / (cnts * r.stat.meanpi0)
        ymin = min(pred.min(), obs.min())
        ymax = max(pred.max(), obs.max())

        plt.plot(bins, obs, color='gray', lw=2)
        plt.plot(bins, pred, color=cols[i], alpha=0.75, lw=1.5)
        i += 1

        plt.xticks(fontsize=20)
        plt.xlabel('cM to nearest {}'.format(focal_dict[focal]), fontsize=22, labelpad=10)
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
        plt.savefig(fig_file, dpi=256)
        # plt.show()
        plt.close()
    # plt.show()


def rsq_plot():
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/rsq_maps/'
    # folder = 'alternateclean'
    # folder = 'partitions'
    folder = 'variants'
    # folder = 'popsclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    d = rdir + folder
    files = ['{}/{}'.format(d, f) for f in os.listdir(d)]
    # files = ['{}/{}'.format(d, f) for f in os.listdir(d)][::-1]

    # labels = ['{}%'.format(p) for p in xrange(1, 11)]
    # labels = [f.split('/')[-1].split('.')[0] for f in files]
    # labels = ['CEU-cons', 'CEU-cons-coding/noncoding', 'LWK-cons', 'LWK-cons-coding/noncoding', 'YRI-cons',
    #           'YRI-cons-coding/noncoding']
    labels = ['YRI-CADD', 'YRI-cons', 'YRI-cons-coding/noncoding']
    arrays = [np.loadtxt(f)[:-2].T for f in files]
    win = np.log10(arrays[0][0])

    plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.09, bottom=0.1, right=1.0, top=1.0, wspace=0.1, hspace=0.12)

    space = 0.013
    shift = -space * 0.5 * len(arrays)
    # shift = 0

    for (i, arr) in enumerate(arrays):
        facecol = clr.colorConverter.to_rgba(cols[i], alpha=0.25)
        edgecol = clr.colorConverter.to_rgba(cols[i], alpha=0.75)
        mdic = dict(markerfacecolor=facecol, markeredgecolor=edgecol, markeredgewidth=1.5, color=edgecol, lw=0.0)
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

    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/variance_explained'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, folder)
    plt.savefig(fig_file, dpi=213)
    # plt.show()


def dfe_plot():
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/final_files/'
    # folder = 'alternateclean'
    # folder = 'popsclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    folder = 'genicpartclean'
    d = rdir + folder
    files = ['{}/{}'.format(d, f) for f in os.listdir(d)]
    # files = ['{}/{}'.format(d, f) for f in os.listdir(d)][::-1]
    rs = [RunStruct(init=f) for f in files]

    # labels = ['{}%'.format(p) for p in xrange(1, 11)]
    # labels = [f.split('/')[-1].split('.')[0] for f in files]
    labels = [['CEU genic cons', 'LWK genic cons', 'YRI genic cons'],
              ['CEU nongenic cons', 'LWK nongenic cons', 'YRI nongenic cons']]

    figwidth = 10
    plt.figure(figsize=(figwidth, 7))
    plt.subplots_adjust(left=0.1, bottom=0.13, right=1, top=1, wspace=0.1, hspace=0.12)

    x = np.arange(rs[0].bsgrid)
    for i in xrange(rs[0].bnum):
        plt.subplot(rs[0].bnum, 1, i+1)
        arr = np.array([r.stat.upmf[i] for r in rs])
        # x = np.arange(len(arr[0]))
        width = 0.8 / len(arr)
        shift = -(0.5 * len(arr) * width)
        for (j, row) in enumerate(arr):
            # plt.plot(x, np.cumsum(row), '-s', color=cols[i])
            plt.bar(left=x + shift, height=row, width=width, align='edge', color=cols[j])
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
    plt.xticks(x, [r'$\mathrm{10^{%.1f}}$' % s for s in np.arange(-4.5, -1.5, 0.5)], fontsize=20)
    plt.xlabel('fitness effects of deleterious mutations', fontsize=22, labelpad=15)

    # leg = plt.legend(labels, prop={'size': 16}, frameon=True, loc='upper left')
    # frame = leg.get_frame()
    # frame.set_facecolor('white')
    # frame.set_alpha(1.0)

    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/dfe'
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, folder)
    # plt.savefig(fig_file, dpi=2560/figwidth)
    plt.show()


def udel_plot():
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/final_files/'
    # folder = 'alternateclean'
    # folder = 'popsclean'
    # folder = 'phyloclean'
    # folder = 'pr90to99clean'
    folder = 'genicpartclean'

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
    plt.savefig(fig_file, dpi=256)


def cons_filters():
    dr = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/cons_filters/'
    ar = [np.loadtxt(dr+f, dtype=str).T for f in os.listdir(dr) if f.startswith('prosimian.')]
    nm = ['{}% threshold'.format(p) for p in xrange(80, 100, 5)]

    width = 0.9 / len(ar)
    shift = -0.4
    x = np.arange(len(ar[0][0])-1)

    plt.figure(figsize=(13, 4))
    plt.subplots_adjust(left=0.06, bottom=0.08, right=1, top=0.97, wspace=0.1, hspace=0.12)
    for (p, a) in zip(nm, ar):
        a = a[1].astype(float)
        plt.bar(left=x+shift, height=a[1:]/a[0], width=width, align='edge', label=p)
        shift += width

    plt.ylabel('fraction remaining', fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0.9, 1.0)

    # plt.xlabel('phylogeny', fontsize=22)
    plt.xticks(x, ar[0][0][1:], fontsize=12)

    plt.legend(ncol=4, loc='upper center', prop={'size': 12})

    plt.show()


def alignment_coverage():
    fname = result_dir + '100way_species_coverage/autosomes.100way.cov.txt'
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

def main():
    # sorted_plot()
    # rsq_plot()
    # dfe_plot()
    # udel_plot()
    # for s in 'syn nonsyn cds_midpoints'.split():
    #     collated_plot(focal=s)
    # collated_plot('syn')
    # cons_filters()
    alignment_coverage()

if __name__ == '__main__':
    main()
