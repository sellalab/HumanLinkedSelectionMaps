from classes.runstruct import ChromStruct, root_dir, human_autosomes
from classes.bkgdcalculator import BkgdMap, BkgdCalculator
from classes.geneticmap import GeneticMap
from classes.annosegments import AnnoSegments
from data_processing.data_tools import calculate_error
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from itertools import izip
import numpy as np
import seaborn

__author__ = 'davidmurphy'


def getidx2dbl():
    INTERP_TAB_FR_MULT = 10.0
    INTERP_TAB_MIN_EXP = -44
    INTERP_TAB_MAX_EXP = 5

    INTERP_TAB_N_ENTRIES = (INTERP_TAB_MAX_EXP - INTERP_TAB_MIN_EXP) * \
                           INTERP_TAB_FR_MULT + 1

    return np.array([idx2dbl(i) for i in xrange(int(INTERP_TAB_N_ENTRIES))])


def idx2dbl(idx):
    """copy of mcvicker's function to convert table indices to floats"""
    INTERP_TAB_FR_MULT = 10.0
    INTERP_TAB_FR_OFFSET = 0.5
    INTERP_TAB_MIN_EXP = -44
    INTERP_TAB_MAX_EXP = 5

    # INTERP_TAB_N_ENTRIES = (INTERP_TAB_MAX_EXP - INTERP_TAB_MIN_EXP) * INTERP_TAB_FR_MULT + 1

    fr = (idx-1) % int(INTERP_TAB_FR_MULT)
    x = (idx-1) / INTERP_TAB_FR_MULT

    if idx < 0:
        raise ValueError('Values below 0 not allowed')
    if idx == 0:
        return 0.0

    return (((INTERP_TAB_FR_OFFSET*fr)/INTERP_TAB_FR_MULT)+INTERP_TAB_FR_OFFSET) * pow(2.0, x+INTERP_TAB_MIN_EXP)


def dbl2idx(dbl):
    """copy of mcvicker's function to convert floats to table indices"""

    if dbl == 0:
        return 0


def cmmb_grad(chrom):
    # use ChromStruct to get gmap file path
    cst = ChromStruct(chrom=chrom)
    # load gmap rates and positions
    pos, cmmb = np.loadtxt(cst.gmap_files, usecols=(0, 1)).T
    # set cmmb==0 positions to 1% of the next lowest cmmb
    cmmb[(cmmb == 0)] = 0.01 * np.unique(cmmb)[1]
    # convert rates to log scale
    log10_cmmb = np.log10(cmmb)
    # find the rage of the log scale of rates, divide into 100 bins
    step_size = (log10_cmmb.max() - log10_cmmb.min()) / 100.0
    bins = np.arange(log10_cmmb.min(), log10_cmmb.max(), step_size)
    # blue channel runs from 100-0%
    red = 0.01 * np.searchsorted(bins, log10_cmmb)
    green = 0.75 * red
    blue = 1 - red
    # create a gradient of colors for each cmmb
    grad = np.column_stack((red, green, blue))

    # safety measure for values just out of bounds
    grad[grad < 0] = 0
    grad[grad > 1] = 1

    return pos, grad


def cons_grad(chrom):
    # use ChromStruct to get cons segments file path
    cst = ChromStruct(chrom=chrom)
    coords = np.load(cst.all_cons_segs)[chrom]
    # tile conserved points on a blank chrom-length zeros array
    carr = np.zeros(shape=cst.chlen, dtype='u1')
    for (i, j) in coords:
        carr[i:j] = 1
    # sum conserved sites every 10kb block across position in the chrom
    block = int(1e4)
    pos = np.arange(0, cst.chlen, block)
    csum = np.array([np.sum(carr[i:i + block], dtype=int) for i in pos])
    # find the range of conserved site sums, divide into 100 bins
    step_size = (csum.max() - csum.min()) / 100.0
    bins = np.arange(csum.min(), csum.max(), step_size)
    # create a blue tined gray scale for conservation density
    green = red = 0.8 * (1 - np.searchsorted(bins, csum) / 100.0)
    blue = 0.25 + (65 * green / 80)
    grad = np.column_stack((red, green, blue))

    return pos, grad


def add_gradient(pos, grad, xmin, xmax, base, height, alpha):
    # mask pos to fit B map range
    msk = (pos >= xmin) & (pos <= xmax)
    pos, grad = pos[msk], grad[msk]
    # define widths as between-position distances for gradient points
    # pp = np.concatenate(([xmin], pos))
    pp = np.concatenate((pos, [xmax]))
    widths = pp[1:] - pp[:-1]
    # widths = np.concatenate((pos[1:] - pos[:-1], [xmax - pos[-1]]))
    # create list of rectangles with appropriate widths for the gradient
    # rects = [Rectangle((x, base), w, height) for (x, w) in zip(pos, widths)]
    rects = [Rectangle((x, base), w, height) for (x, w) in zip(pp[:-1], widths)]
    # create a "collection" from the rectangles
    collection = PatchCollection(rects)
    # define features of grad rectangles
    collection.set_alpha(alpha)
    collection.set_facecolor(grad)
    collection.set_linewidth(0)

    return collection


def compare_dists(ti, an, sc):
    cst = ChromStruct(chrom='chr1')

    vals = []
    for ch in cst.chroms:
        bxct = '{}/precalc/b_exact/{}.b.exact.rand.{}.{:.6f}.npy'
        sx, sy = np.load(bxct.format(root_dir, ch, an, ti)).T

        bdir = 'bigB{}'.format(int(sc))
        cst.bdir = bdir
        midx = np.where(cst.bdfe[0] == ti)[0][0]
        bbx, bby = BkgdMap(cst.bmap_files[an][midx], False).bmap
        # bbx, bby = load_bmap(cst.bmap_files[an][midx], sc, False).T
        bi = np.minimum(np.searchsorted(bbx, sx), bbx.size - 1)

        bdir = 'littleB{}'.format(int(sc))
        cst.bdir = bdir
        lbx, lby = BkgdMap(cst.bmap_files[an][midx]).bmap
        # lbx, lby = load_bmap(cst.bmap_files[an][midx], sc).T
        li = np.minimum(np.searchsorted(lbx, sx), lbx.size - 1)

        vals.append(np.column_stack((sy, bby[bi], lby[li])))

    sy, bby, lby = np.concatenate(vals).T

    plt.figure(figsize=(13, 7))
    coeff = r' $t=10^{%.1f}$' % np.log10(ti)
    title = 'B distributions: {}'.format(coeff)
    plt.title(title, fontsize=16)
    traits = dict(lw=0.8, cumulative=1, histtype='step', bins=sc, normed=1)
    plt.hist(sy, label=r'$B_{exact}$', color='darkturquoise', **traits)
    plt.hist(bby, label=r'$B_{mcvicker}$', color='fuchsia', **traits)
    plt.hist(lby, label=r'$B_{new}$', color='darkorange', **traits)
    plt.xticks(fontsize=16)
    plt.xlabel('B values', fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel('proportion', fontsize=16)
    plt.legend(prop={'size': 14}, loc='upper left')
    plt.show()


def compare_error(ti, an, sc, rel_err=True):
    # cutoff threshold for outliers
    threshold = 2.0 / sc  # (x2 for two steps where error is invoked)
    bounds = 2 * threshold
    nbins = 100.0
    # divide errors into 100 bins inside the outlier range
    bins = np.arange(-bounds, bounds, bounds/nbins)

    # get exactB, littleB, bigB and errors for each map across chroms
    # try:
    #     err = [calc_err(ch, ti, an, sc, rel_err)[-1] for ch in human_autosomes]
    # except ValueError:
    #     print 'print WHAT?!'
    # aut = [human_autosomes[0]] + list(human_autosomes[2:])
    # aut = human_autosomes
    aut = human_autosomes[6:]
    # # DMEL MODE
    # aut = dmel_autosomes
    err = [calc_err(ch, ti, an, sc, rel_err)[-1] for ch in aut]
    err = np.concatenate(map(np.column_stack, err)).T

    # get the fraction of errors past a range of benchmarks
    ibar = []
    for frac in threshold * 2**np.arange(5):
        pct = [100 * np.sum(abs(r) > frac) / float(r.size) for r in err]
        ibar.append(pct)
    ibar = np.array(ibar)
    ibar[:-1] -= ibar[1:]

    # use the same plot bounds for each set of values
    ibnd = [abs(r) <= bounds for r in err]
    hist = [np.histogram(r[i], bins) for r, i in izip(err, ibnd)]

    # plot map comparisons and annotations
    fig = plt.figure(figsize=(6.4, 4.8))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.88, top=0.94,
                        wspace=0.15)
    ax1 = fig.add_subplot(1, 5, (1, 4))
    # ax1 = fig.add_subplot(111)
    # map_names = 'littleB bigB'.split()
    # # DMEL MODE
    # map_names = ['dmel100']
    map_names = ['littleB']
    map_cols = ['fuchsia', 'darkorange']
    line_objs = []
    sdict = dict(linewidth=0.8, alpha=0.75)
    for (i, hs) in enumerate(hist):
        # if i > 0:
            # break
            # pass
        freq = hs[0].astype(float) / hs[0].sum()
        bins = hs[1][1:]
        line, = plt.step(bins, freq, color=map_cols[i], **sdict)
        line_objs.append(line)

    # draw vertical lines at expected max error bounds
    plt.axvline(-threshold, color='k', ls='--', alpha=0.5, lw=0.6)
    xtol = plt.axvline(threshold, color='k', ls='--', alpha=0.5, lw=0.6)
    map_names.append('bound')
    line_objs.append(xtol)

    # create the plot title
    pex, tex = (1.0 / sc), np.log10(ti)
    rcalc = 'relative\ error' if rel_err else 'additive\ error'
    precs = r'$\mathrm{\epsilon=%.3f}$' % pex
    coeff = r' $\mathrm{(t=10^{%.1f})}$' % tex
    title = precs + coeff
    plt.suptitle(title, fontsize=11)

    # format axes
    x_points = [x*threshold for x in -2, -1, 0, 1, 2]
    x_labels = [r'$\mathrm{%d\epsilon}$' % x if x != 0 else '0'
                for x in -4, -2, 0, 2, 4]
    plt.xticks(x_points, x_labels, fontsize=11)
    plt.xlabel(r'$\mathrm{%s}$' % rcalc, fontsize=11)
    plt.yticks(fontsize=11)
    plt.ylabel('fraction', fontsize=12)
    ax1.legend(line_objs, map_names, loc=2, fontsize=10)

    # outlier bar plots
    ax2 = fig.add_subplot(155)
    # bar_objs = []
    fracs = []
    bots = [0, 0]
    xs = [-0.2, 0.8]
    cols = np.array(['fuchsia blueviolet darkviolet purple indigo'.split(),
                     'darkorange tomato crimson firebrick darkred'.split()])
    for (i, pct) in enumerate(ibar):
        for (j, pc) in enumerate(pct):
            # plot bar segment for outliers if segment >0.01%
            if pc > 0.01:
                plt.bar(j, pc, bottom=bots[j], color=cols[j, i])
                txt = '>{:.0f}'.format(2**(i + 1)) + r'$\mathrm{\epsilon}$'
                plt.text(x=xs[j], y=bots[j] + 0.5 * pc, s=txt, fontsize=10)

                # add legend object if its the first pass
                if i == 0:
                    fracs.append(ibar[:, j].sum())

                # reset bottom coordinates for bars
                bots[j] += pc

        # bar1, = plt.bar(0, pct[0], bottom=bot1, color=cols[0, i])
        # bar2, = plt.bar(1, pct[1], bottom=bot2, color=cols[1, i])
            # if pct[1] > 0.01:
            #     plt.text(x=0.8, y=bot2 + 0.5 * pct[1], s=txt, fontsize=10)
        # if i == 0:
        #     bar_objs = [bar1, bar2]
        #     fracs = [ibar[:, 0].sum(), ibar[:, 1].sum()]

        # bot1 += pct[0]
        # bot2 += pct[1]

    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    # plt.title(r'$\mathrm{err>2\epsilon}$', fontsize=10)
    plt.xticks([0, 1], ['{:.2g}%'.format(x) for x in fracs], fontsize=11)
    plt.yticks(fontsize=11)
    plt.ylabel(r'$\mathrm{err>2\epsilon\ (\%)}$', fontsize=10)

    # bar_objs.append(xtol)
    # ax1.legend(bar_objs, map_names, loc=2, fontsize=10)
    plt.grid('off')
    # plt.show()

    rcalc = 'rel' if rel_err else 'add'
    tstr = '{}_t=1e{:.1f}_eps=1e{:.1f}_{}_.png'
    # # DMEL MODE
    # tstr = '{}_t=1e{:.1f}_eps=1e{:.1f}_{}_dmel.png'
    token = tstr.format(rcalc, tex, np.log10(pex), an)
    folder = 'murphy_et_al/figures/precision_tests/err_dist'
    fname = '{}/{}/{}'.format(root_dir, folder, token)
    plt.savefig(fname, dpi=256)
    plt.close()


def chrom_error(ti, an, sc, rel_err=True):
    # cutoff threshold for outliers
    threshold = 2.0 / sc  # (x2 for two steps where error is invoked)
    bounds = 2.5 / sc
    nbins = 100.0
    # divide errors into 100 bins inside the outlier range
    bins = np.arange(-bounds, bounds, bounds/nbins)

    # get exactB, littleB, bigB and errors for each map across chroms
    err = [calc_err(ch, ti, an, sc, rel_err)[-1] for ch in human_autosomes]
    # err = np.concatenate(map(np.column_stack, err)).T

    # get the fraction of errors outside of expected bounds
    lfrac = [np.sum(abs(r[0]) > threshold) / float(r[0].size) for r in err]
    bfrac = [np.sum(abs(r[1]) > threshold) / float(r[1].size) for r in err]
    frac = [lfrac, bfrac]

    # plot error/chrom if there is any >2epsilon error on any chrom
    if any(any(fr) for fr in frac):
        plt.figure(figsize=(9, 5))
        plt.subplots_adjust(left=0.08, bottom=0.16, right=1, top=0.94)
        map_names = 'littleB bigB'.split()
        map_cols = ['fuchsia', 'darkorange']
        sdict = dict(linewidth=0.8, alpha=0.75)
        width = 0.4
        xvals = np.arange(22)

        for (i, fr) in enumerate(frac):
            if any(fr):
                plt.bar(left=xvals-width*i, height=fr, width=width,
                        color=map_cols[i], label=map_names[i], **sdict)

        plt.xticks(xvals, human_autosomes, rotation=90)
        plt.xlabel('chromosome')
        plt.ylabel(r'$\mathrm{fraction\ >2\epsilon}$')
        plt.legend()

        # create the plot title
        pex, tex = (1.0 / sc), np.log10(ti)
        rcalc = 'relative\ error' if rel_err else 'additive\ error'
        precs = r'$\mathrm{%s:\ \epsilon=%.3f}$' % (rcalc, pex)
        coeff = r' $\mathrm{(t=10^{%.1f})}$' % tex
        title = precs + coeff
        plt.title(title)

    # rcalc = 'rel' if rel_err else 'add'
    # token = '{}_t=1e{:.1f}_eps=1e{:.1f}.png'.format(rcalc, tex, np.log10(pex))
    # folder = 'murphy_et_al/figures/precision_tests/err_dist'
    # fname = '{}/{}/{}'.format(root_dir, folder, token)
    # plt.savefig(fname, dpi=256)
    # plt.close()


def compare_maps(ch, ti, an, sc, rel_err=True):
    threshold = 2.0 / sc

    # get exactB, estimates from both maps and errors for each map
    b_xct, b_est, b_err = calc_err(ch, ti, an, sc, rel_err)

    # get outlier sites where the error is over 5 scale units off
    iout = [abs(r) > threshold for r in b_err]
    outliers = [(b_xct[0][i], b[i]) for (b, i) in izip(b_est, iout)]

    # if there are any outliers in either estimated B map, show on plot
    if any(r[0].any() for r in outliers):

        # plot map comparisons and annotations
        fig = plt.figure(figsize=(13, 7))
        ax = fig.add_subplot(4, 1, (1, 3))
        plt.subplots_adjust(left=0.08, bottom=0.1, right=1, top=0.94)

        # plot B exact
        legend_labs = []
        legend_objs = []
        sx, sy = b_xct
        bxt, = plt.plot(sx, sy, lw=0.8, color='k')
        legend_objs.append(bxt)
        legend_labs.append('exactB')

        # plot each B estimate
        map_cols = ['fuchsia', 'darkturquoise']
        # map_names = 'littleB bigB'.split()
        # DMEL MODE
        map_names = ['dmel100']
        for (i, bs) in enumerate(b_est):
            bplot, = plt.plot(sx, bs, lw=0.8, color=map_cols[i])
            legend_objs.append(bplot)
            legend_labs.append(map_names[i])

        # plot red circles around regions where B estimate is far off B exact
        err_cols = ['red', 'lime']
        err_mark = ['s', 'o']
        rdict = dict(facecolors='none', linewidth=1, alpha=0.75)
        for (i, er) in enumerate(outliers):
            errplot = plt.scatter(er[0], er[1], edgecolors=err_cols[i],
                                  marker=err_mark[i], **rdict)
            legend_objs.append(errplot)
            legend_labs.append(map_names[i] + ' outlier')

        # plt.show()

        alph = 1
        # create blue-to-orange heat map of recombination rates
        pos, grad = cmmb_grad(chrom=ch)
        col = add_gradient(pos, grad, sx.min(), sx.max(), -0.2, 0.2, alph)
        ax.add_collection(col)

        # create green-to-red heat map of conservation density
        pos, grad = cons_grad(chrom=ch)
        col = add_gradient(pos, grad, sx.min(), sx.max(), -0.41, 0.2, alph)
        ax.add_collection(col)

        # label gradients
        gdict = dict(marker='s', markersize=20, ls='none', alpha=alph)
        # double-sided block label for cmmb
        m1, = plt.plot([], [], c=(0, 0, 1), fillstyle='left', **gdict)
        m2, = plt.plot([], [], c=(1, 0.75, 0), fillstyle='right', **gdict)
        cmmb = (m1, m2)
        # double-sided block label for cons density
        m3, = plt.plot([], [], c=(0.8, 0.8, 0.9), fillstyle='left', **gdict)
        m4, = plt.plot([], [], c=(0, 0, 0.25), fillstyle='right', **gdict)
        cons = (m3, m4)
        spacer, = plt.plot([], [], alpha=0)
        # create legend
        legend_objs += [spacer, cons, cmmb]
        legend_labs += ['', 'conservation level', 'log10(cM/Mb)']

        ldict = dict(numpoints=1, loc=2, fontsize=14, ncol=3, labelspacing=1.1)
        plt.legend(legend_objs, legend_labs, **ldict)
        # leg = plt.legend(legend_objs, legend_labs, **ldict)
        # frame = leg.get_frame()
        # frame.set_facecolor('white')
        #
        rcalc = 'relative\ error' if rel_err else 'additive\ error'
        coeff = r' $\mathrm{%s\ (t=10^{%.1f})}$' % (rcalc, np.log10(ti))
        title = ch + r': $\mathrm{B_{exact}\ and\ B_{estimate}}$ ' + coeff
        plt.title(title, fontsize=16)
        plt.xticks(fontsize=16)
        plt.xlabel('chromosome position (bp)', fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylabel('B values', fontsize=16)
        plt.ylim(-0.45, 1.25)
        # plt.ylim(0.6, 1.1)
        plt.show()

    return None


def calc_err(ch, ti, an, sc, rel_err=True):
    # load B exact

    # # new B_exact from interpolator endpoints
    # mdir = root_dir + '/precalc/exactPyB_vs_exactMcvB'
    # fname = mdir + '/{}.PyVsC.0.000032.npy'
    # b_xct = np.load(fname.format(ch))[:, :2].T

    # # DMEL MODE
    # fexact = '{}/precalc/dmel_b_exact/{}.b.exact.rand.{}.{:.6f}.npy'

    fexact = '{}/precalc/b_exact_1e5/{}.b.exact.rand.{}.{:.6f}.npy'
    b_xct = np.load(fexact.format(root_dir, ch, an, ti)).T

    # interpolate B values in each map at B exact x-coords
    b_err = []
    b_est = []

    # # --> FOR DMEL ONLY!
    # coefs = 10 ** np.arange(-1.5, -6, -1)
    # cst = ChromStruct(chrom='2L',
    #                   gmap='comeron', mask='', ncon='', cons='',
    #                   outg='', neut='dmel', phase='', token='MISSING',
    #                   focal='exonicNS', rm_annos=(), bkgd_scale=100.0,
    #                   bs_annos=('exon', 'longintron', 'intergenic', 'UTR'),
    #                   bs_dfe=tuple([coefs] * 4),
    #                   bmap_dir='dmel',
    #                   cs_annos=('exonicNS', 'intronic', 'UTR'),
    #                   cs_dfe=tuple([coefs] * 3),
    #                   cmap_dir='dmel', chroms=dmel_autosomes,
    #                   methods=('Nelder-Mead', 'COBYLA', 'Powell'))
    # cst.fixed.u_fix = 3.5e-9
    # # <-- FOR DMEL ONLY!

    # map_names = 'littleB bigB'.split()
    map_names = ['littleB']
    # # DMEL MODE
    # map_names = ['dmel']
    for bm in map_names:
        # load mcvicker bmaps and compare with B exact
        cst = ChromStruct(ch, bdir='{}{}'.format(bm, int(sc)))
        midx = np.where(cst.bdfe[0] == ti)[0][0]
        bexp = bool('littleB' in bm)
        bfile = cst.bmap_files[an][midx]
        bmp = BkgdMap(bfile, bexp)

        # interpolate b values at the exact bmap sample sites
        # b_int = np.interp(b_xct[0], bmp.pts, bmp.b)
        b_interp = bmp.interp_b(b_xct[0])
        # calculate error between exact and estimated b
        err = bmap_err(b_xct[1], b_interp, rel_err)

        b_est.append(b_interp)
        b_err.append(err)

    return b_xct, b_est, b_err


def bmap_err(b_exact, b_estimate, relative=True):
    if relative:
        return (b_exact - b_estimate) / b_exact
    else:
        return b_exact - b_estimate


def calc_error(fexact, fbkgd, rel_err=True, little_b=True):
    # load exact B values
    b_xct = np.load(fexact).T
    # load bkgd map
    bmp = BkgdMap(fbkgd, little_b)
    # get bkgd map values at exact B positions
    b_interp = bmp.interp_b(b_xct[0])

    # calculate relative/additive error between exact and map B values
    err = bmap_err(b_xct[1], b_interp, rel_err)

    return np.column_stack((b_xct.T, b_interp, err))


def hist_error(eps, nbins, rlist):
    tol = 2.0 * eps  # tolerance range = 2 x epsilon
    bnd = 2.0 * tol  # plot bounds = 4 x epsilon
    inc = 2.0 * bnd / nbins  # bin increments according to number of bins
    bns = np.arange(-bnd, bnd+inc, inc)

    hlist = []  # list for histogram counts of error values
    olist = []  # list for outliers
    out = tol * 2.0 ** np.arange(5)
    for r in rlist:
        # make a histogram of errors within +/-4 x epsilon
        histo = np.histogram(r, bns)
        hlist.append(histo)
        # count outliers above log range of thresholds
        outlr = np.array([100.0 * np.sum(abs(r) > o) / r.size for o in out])
        outlr[:-1] -= outlr[1:]
        olist.append(outlr)

    return hlist, olist


def plot_error_dist(t, eps, hlist, olist, names, token):
    tol = 2.0 * eps  # tolerance range = 2 x epsilon
    map_cols = ['fuchsia', 'darkorange', 'darkturquoise']
    bar_cols = np.array([
        ['fuchsia', 'blueviolet', 'darkviolet', 'purple', 'indigo'],
        ['darkorange', 'tomato', 'crimson', 'firebrick', 'darkred']])
    sdict = dict(linewidth=0.8, alpha=0.75)

    # plot error distribution
    fig = plt.figure(figsize=(6.4, 4.8))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.88, top=0.94,
                        wspace=0.15)
    # format title
    pex, tex = eps, np.log10(t)
    # rcalc = 'relative\ error'
#    tk = r'$\mathrm{%s}$: ' % token
    tk = token

    precs = r'$\mathrm{\epsilon=%.1g}$' % pex
    # coeff = r' $\mathrm{(t=10^{%.1f})}$' % tex
    coeff = r' $\mathrm{(s=10^{%.1f})}$' % tex
    title = tk + precs + coeff
    plt.suptitle(title, fontsize=11)

    # plot distribution of errors in leftmost 80% of plot area
    line_objs = []  # store line handles
    ax1 = fig.add_subplot(1, 5, (1, 4))
    for (i, h) in enumerate(hlist):
        freq = h[0].astype(float) / h[0].sum()
        bins = h[1][1:]
        line, = ax1.step(bins, freq, color=map_cols[i], **sdict)
        line_objs.append(line)
    ax1.axvline(-tol, color='k', ls='--', alpha=0.5, lw=0.6)
    xtol = ax1.axvline(tol, color='k', ls='--', alpha=0.5, lw=0.6)
    names.append(r'$\mathrm{tolerance}$')
    line_objs.append(xtol)

    # plot x labels in units of epsilon
    x_points = [x*tol for x in -2, -1, 0, 1, 2]
    x_labels = [r'$\mathrm{%d\epsilon}$' % x if x != 0 else r'$\mathrm{0}$'
                for x in -4, -2, 0, 2, 4]
    plt.xticks(x_points, x_labels, fontsize=12)
    plt.xlabel('relative error', fontsize=12)
    plt.yticks(fontsize=11)
    plt.ylabel('fraction', fontsize=12)
    ax1.legend(line_objs, names, loc=2, fontsize=10)

    # plot bar plots of outlier values stratified into epsilon * 2^x bins
    ax2 = fig.add_subplot(155)
    fracs = []
    bots = np.zeros(len(olist))  # bar plot bottoms
    xs = np.arange(len(olist)) + 0.4  # bar segment text label positions
    for (i, o) in enumerate(olist):
        fracs.append(o.sum())
        for (j, p) in enumerate(o):
            # plot bar segment for outliers if segment >0.01%
            if p > 0.01:
                # NOTE: i shifted the bars from i+0.1 to i+0.5 because
                plt.bar(i+0.5, p, bottom=bots[i], color=bar_cols[i, j])
                txt = '>{:.0f}'.format(2**(j + 1)) + r'$\mathrm{\epsilon}$'
                xi, yi = xs[i], bots[i] + 0.5 * p
                plt.text(x=xi, y=yi, s=txt, fontsize=12, color='white')
                bots[i] += p

    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    # plt.title(r'$\mathrm{err>2\epsilon}$', fontsize=10)
    plt.xticks([0.5, ], ['{:.2g}%'.format(x) for x in fracs], fontsize=11)
    plt.yticks(fontsize=11)
    plt.xlim(0, 1)
    plt.ylabel(r'$\mathrm{err>2\epsilon\ (\%)}$', fontsize=12)
    plt.grid('off')
    # plt.show()


def main_1():
    a = 'primate_cons95_Segments'
    t = 10**-4.5
    eps = 0.02
    cst = ChromStruct('chr1')
    x_str = root_dir + '/precalc/b_exact_1e5/{}.b.exact.rand.{}.{}.npy'
    bs1 = root_dir + '/precalc/littleB50/AA_Map_10kb.{}.{}.t{}.bkgd'
    bs2 = root_dir + '/precalc/littleB50-old/AA_Map_10kb.{}.{}.t{}.bkgd'

    chroms = cst.chroms[10:]
    t_str = '{:.6f}'.format(t)
    ex = [x_str.format(c, a, t_str) for c in chroms]
    t_str = '{:.8f}'.format(t)
    b1 = [bs1.format(c, a, t_str) for c in chroms]
    b2 = [bs2.format(c, a, t_str) for c in chroms]

    r1 = np.concatenate([
        calc_error(fx, fb)[:, -1] for (fx, fb) in izip(ex, b1)])
    r2 = np.concatenate([
        calc_error(fx, fb, 50.0)[:, -1] for (fx, fb) in izip(ex, b2)])

    hl, ol = hist_error(eps, 200, [r1, r2])

    plot_error_dist(t, eps, hl, ol, ['bkgd-new', 'bkgd-old'])


def main_2(folder, tok, ch, samples):
    an = 'primate_cons95_Segments'
    t = 10 ** -4.5
    eps = 0.02
    tol = 2.0 * eps

    # bf = root_dir + '/precalc/littleB50Debug_2/{}'.format(folder)
    # bfnam = '{}/AA_Map_10kb.{}.{}.t{:.8f}.exactB.txt'.format(bf, ch, an, t)
    bf = root_dir + '/precalc/littleB50Debug'
    ofnam = '{}/{}.exactB.calculated.t{:.8f}.npy'.format(bf, ch, t)
    bx, by, smp = np.load(ofnam).T

    # bx, by = np.loadtxt(bfnam).T
    # by = np.exp(by)
    # i = np.unique(np.random.randint(0, len(bx), samples))
    # bx, by = bx[i], by[i]
    # bfnam = '{}/AA_Map_10kb.{}.{}.t{:.8f}.bkgd'.format(bf, ch, an, t)
    # bmp = BkgdMap(bfnam)
    # bx = bmp.mids
    # by = bmp.interp_b(bx)

    # cst = ChromStruct(ch)
    # gfnam = cst.gmap_files
    # gmp = GeneticMap(ch, gfnam, 0.01)

    # afnam = '{r}/data/bsanno/{a}/{f}/{c}.{a}.bed'.format(
    #     r=root_dir, a=an, f=folder, c=ch)
    # # afnam = '{r}/data/bsanno/{a}/{c}.{a}.bed'.format(r=root_dir, a=an, c=ch)
    # ano = AnnoSegments(ch, gmp, afnam)
    # bkc = BkgdCalculator(t, cst.fixed.u_fix, ch, gmp, ano)
    # smp = np.exp(np.array([bkc.bkgd_sum_chrom(y=y) for y in bx]))
    r1 = bmap_err(smp, by, 1)  # use standard error plot

    # EXACT VS ESTIMATE ALONG CHROMOSOME PLOT
    sdict = dict(linewidth=0.8, alpha=0.75)
    token = r'$\mathrm{%s}:\ $' % tok
    pdir = 'murphy_et_al/figures/precision_tests/modify_code_and_data'
    plt.figure(figsize=(12, 4.8))
    plt.subplots_adjust(left=0.08, bottom=0.12, right=0.97, top=0.94)
    pex, tex = eps, np.log10(t)
    precs = r'$\mathrm{\epsilon=%.3f}$' % pex
    coeff = r' $\mathrm{(t=10^{%.1f})}$' % tex
    plt.suptitle(token + precs + coeff, fontsize=11)
    plt.plot(1e-7 * bx, by, label='B_estimate', color='fuchsia', **sdict)
    plt.plot(1e-7 * bx, smp, label='B_exact', color='darkorange', **sdict)
    plt.xlabel(r'$\mathrm{position\ (x10Mbp)}$', fontsize=11)
    plt.xticks(fontsize=11)
    plt.ylabel(r'$\mathrm{B\ values}$', fontsize=11)
    plt.yticks(fontsize=11)
    plt.legend()
    # fimg = '{}/{}/{}.{}.bmap.png'.format(root_dir, pdir, ch, folder)
    # plt.savefig(fimg, dpi=256)
    # plt.close()

    # ERROR MEASUREMENTS PLOTTED ALONG CHROMOSOME
    fig = plt.figure(figsize=(12, 4.8))
    ax1 = fig.add_subplot(5, 1, (1, 4))
    plt.subplots_adjust(left=0.08, bottom=0.12, right=0.97, top=0.94, hspace=0)
    pex, tex = eps, np.log10(t)
    precs = r'$\mathrm{\epsilon=%.3f}$' % pex
    coeff = r' $\mathrm{(t=10^{%.1f})}$' % tex
    plt.suptitle(token + precs + coeff, fontsize=11)
    sct = plt.scatter(1e-7 * bx, r1, facecolors='none', edgecolors='dodgerblue',
                      **sdict)
    y_rng = list(-2 ** np.arange(4, 0, -1)) + [0] + list(2 ** np.arange(1, 5))
    y_points = [tol * y for y in y_rng]
    y_labels = [r'$\mathrm{%d\epsilon}$' % y if y != 0 else '0' for y in y_rng]
    plt.yticks(y_points, y_labels, fontsize=10)
    plt.ylabel(r'$\mathrm{relative\ error}$', fontsize=11)
    plt.ylim(y_points[0], y_points[-1])
    plt.xlabel(r'$\mathrm{position\ (x10Mbp)}$', fontsize=11)
    plt.xticks(fontsize=11)
    # plt.legend()

    # GRADIENT PLOTS FOR CMMB AND CONS
    ax2 = fig.add_subplot(515)
    # create green-to-red heat map of conservation density
    pos, grad = cmmb_grad(chrom=ch)
    col = add_gradient(pos, grad, bx.min(), bx.max(), 0, 0.5, 1.0)
    ax2.add_collection(col)
    # create green-to-red heat map of conservation density
    pos, grad = cons_grad(chrom=ch)
    col = add_gradient(pos, grad, bx.min(), bx.max(), 0.5, 0.5, 1.0)
    ax2.add_collection(col)
    ax2.set_xlim(bx.min(), bx.max())
    ax2.set_ylim(0.0, 1.0)
    # label gradients
    gdict = dict(marker='s', markersize=20, ls='none', alpha=1.0)
    # double-sided block label for cmmb
    m1, = plt.plot([], [], c=(0, 0, 1), fillstyle='left', **gdict)
    m2, = plt.plot([], [], c=(1, 0.75, 0), fillstyle='right', **gdict)
    cmmb = (m1, m2)
    # double-sided block label for cons density
    m3, = plt.plot([], [], c=(0.8, 0.8, 0.9), fillstyle='left', **gdict)
    m4, = plt.plot([], [], c=(0, 0, 0.25), fillstyle='right', **gdict)
    cons = (m3, m4)
    # create legend
    legend_objs = [sct, cons, cmmb]
    legend_labs = ['relative error', 'conservation level', 'log10(cM/Mb)']
    ax1.legend(legend_objs, legend_labs, ncol=3)
    plt.yticks([])
    plt.xlabel(r'$\mathrm{position\ (x10Mbp)}$', fontsize=11)
    plt.xticks(fontsize=11)
    # fimg = '{}/{}/{}.{}.errmap.png'.format(root_dir, pdir, ch, folder)
    # plt.savefig(fimg, dpi=256)
    # plt.close()

    # ERROR DISTRIBUTION PLOT
    hl, ol = hist_error(eps, 200, [r1])
    plot_error_dist(t, eps, hl, ol, ['bkgd-single'], token)
    # fimg = '{}/{}/{}.{}.errdist.png'.format(root_dir, pdir, ch, folder)
    # plt.savefig(fimg, dpi=256)
    # plt.close()

    plt.show()


def main_3():
    t = 10**-4.5
    u = 7.4e-08
    ch = 'chr1'
    ln10 = np.log(10)
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(ch, gmp, coords=[[10, 20]])
    bkc = BkgdCalculator(t, u, ch, gmp, ano)

    # blocks are lognormal with mu=1.4
    # inter-block gaps are lognormal with mu=2.8
    # recombination rates on blocks are lognormal with mu=-8.7 M/bp

    trials = 1000
    samples = 1000
    point = 0.0

    # random sample of block sizes (blocks sizes are approx log-normal)
    blk_mu = 2.0
    blks = np.random.lognormal(ln10 * blk_mu, ln10 * 0.5, samples).astype(int)
    blks = np.maximum(1, blks)

    # random gaps go in between blocks
    gap_mu = 2.0
    gaps = np.random.lognormal(ln10 * gap_mu, ln10 * 0.5, samples).astype(int)

    # random sample of recombination rates (rates are also approx log-normal)
    b_r = np.random.lognormal(ln10 * -8.7, ln10 * 0.6, samples)
    # gap rates are not as low as block rates on average
    g_r = np.random.lognormal(ln10 * -7.9, ln10 * 0.6, samples)

    # interlace blk and gap segments
    segments = np.vstack((gaps, blks)).flatten('F')
    rates = np.vstack((g_r, b_r)).flatten('F')
    r_s = segments*rates  # convert all segments to genetic map units
    r_b = np.cumsum(r_s)[1::2]  # get the ends of conserved block segments
    r_a = r_b - r_s[1::2]  # get the starts of conserved block segments

    # data for each of the b-calculation types (sum and integral)
    rsites = [s + np.arange(b) * r for (s, b, r) in izip(r_a, blks, b_r)]
    rblks = [(a, b, r) for (a, b, r) in izip(r_a, r_b, b_r)]

    log_space = np.logspace(-8, -2, 1000)
    dist = np.concatenate((-log_space, log_space))
    # dist = log_space
    d = point + np.random.choice(dist, trials, False)
    d.sort()
    bs = np.zeros(shape=(2, trials))
    for (i, g) in enumerate(d):
        b_sum = sum(bkc.bkgd_sum(point, g+rs) for rs in rsites)
        b_int = sum(bkc.bkgd_intgr(point, g+a, g+b, r) for (a, b, r) in rblks)
        bs[0, i] = b_sum
        bs[1, i] = b_int

    bs = np.exp(bs)

    plt.figure(figsize=(13, 7))
    # xp = np.floor(np.log10(brate))
    # rmdr = np.log10(brate) - xp
    # coef = 10**rmdr
    # title = r'$\mathrm{block\ length=%dbp}$' + '\n' + \
    #         r'$\mathrm{\ rec.\ rate=%.2f\cdot10^{%.1f}\ M/bp}$'
    # plt.title(title % (blen, coef, xp))

    # bx = np.log10(abs(d-point))
    bx = d
    plt.plot(bx, bs[0], marker='s', label='sum method', alpha=0.75,
             color='darkorange')
    plt.plot(bx, bs[1], marker='o', label='integral method', alpha=0.75,
             color='purple')
    plt.xlabel('(log10) Morgans from block', fontsize=11)
    plt.xticks(fontsize=11)
    plt.ylabel('B values', fontsize=11)
    plt.yticks(fontsize=11)
    # plt.scatter(bx, bs[0], marker='s', label='sum method',
    #             edgecolors='dodgerblue', lw=0.8)
    # plt.scatter(bx, bs[1], marker='o', label='integral approximation',
    #             facecolors='none', edgecolors='darkorange', lw=0.8)
    plt.legend()

    plt.figure(figsize=(13, 7))
    err = bmap_err(bs[0], bs[1], False)
    plt.scatter(bx, err, label='B value error', facecolors='none',
                edgecolors='dodgerblue', lw=1.0)
    plt.legend()

    plt.show()


def main_4():
    t = 10 ** -4.5
    u = 7.4e-08
    ch = 'chr1'
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(ch, gmp, coords=[[10, 20]])
    bkc = BkgdCalculator(t, u, ch, gmp, ano)
    n = 150
    pt = 0.0  # neutral position

    for b_len in 5*10.0**np.arange(4):
        for b_rate in 10.0**np.arange(-12, -2, 2):
            b_pos = np.arange(b_len) * b_rate
            bm1, bm2 = (0.0, b_len * b_rate)
            dist = np.logspace(-15, 0, n)
            bi = np.zeros(shape=(dist.size, 2))
            for (i, d) in enumerate(dist):
                c = pt+d
                bsum = bkc.bkgd_sum(pt, c+b_pos)
                bint = bkc.bkgd_intgr(pt, c+bm1, c+bm2) * b_rate
                bi[i] = bsum, bint

            bs = np.exp(bi.T)
            bx = np.log10(dist)
            # xp = np.floor(np.log10(brate))
            # rmdr = np.log10(brate) - xp
            # coef = 10**rmdr
            plt.figure(figsize=(6.4, 4.8))
            title = r'$\mathrm{block\ length=%dbp}$' + '\n' + \
                    r'$\mathrm{\ rec.\ rate=10^{%.1f}\ M/bp}$'
            plt.title(title % (b_len, np.log10(b_rate)))
            plt.plot(bx, bs[0], marker='s', label='sum method', alpha=0.75,
                     color='darkorange', lw=0.75, ms=5)
            plt.plot(bx, bs[1], marker='o', label='integral method', alpha=0.75,
                     color='purple', lw=0.75, ms=5)
            plt.xlabel('(log10) Morgans from block', fontsize=11)
            plt.xticks(fontsize=11)
            plt.ylabel('B values', fontsize=11)
            plt.yticks(fontsize=11)
            plt.legend()

        plt.show()


def main_5():
    t = 10 ** -4.5
    u = 7.4e-08
    ch = 'chr1'
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(ch, gmp, coords=[[10, 20]])
    # ano = AnnoSegments(ch, gmp, afile=cst.bs_targets.values()[0])
    bkc = BkgdCalculator(t, u, ch, gmp, ano)

    tf = root_dir + '/precalc/littleB50Debug_2/blk.intg.lookup.table.txt'
    with open(tf, 'r') as f:
        f.readline()
        ar = []
        for line in f:
            line = line.split()
            x = float(line[0].split(',')[0].split('=')[1])

            al = []
            for li in line[1:]:
                li = li.split(',')
                # block length to test
                b_len = 25.0
                # y is the distance to the end of segment in M
                y = float(li[0].split('=')[1])
                # assume the segment is 2bp long and calculate M/Mb:
                m_mb = y / b_len
                # correct z values by multiplying by M/Mb and dividing by udel
                z = float(li[1].split('=')[1]) * (-u / m_mb)
                # calculated z value using integral method
                z_calc_int = bkc.bkgd_intgr(0.0, x, x+y) / m_mb
                # calculate z value using sum method assuming 2 bases
                z_calc_sum = bkc.bkgd_sum(0.0, x + np.arange(b_len)*m_mb)
                al.append((x, y, z, z_calc_int, z_calc_sum))

            ar.append(al)
    ar = np.array(ar)

    td = '/Users/davidmurphy/GoogleDrive/linked_selection/lsm_python/temp_data'
    farray = td + '/lkup.tbl.npy'
    np.save(farray, ar)
    ar = np.load(farray)

    a = ar[0]
    plt.plot(a[:, 1], np.exp(a[:, 2]), label='table')
    plt.plot(a[:, 1], np.exp(a[:, 4]), label='exact')
    plt.xscale('log')
    plt.legend()
    plt.show()


def main_6():
    """
    ======================
    3D surface (color map)
    ======================

    Demonstrates plotting a 3D surface colored with the coolwarm color map.
    The surface is made opaque by using antialiased=False.

    Also demonstrates using the LinearLocator and custom formatting for the
    z axis tick labels.
    """

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    # # Make data.
    # x = np.arange(-5, 5, 0.25)
    # y = np.arange(-5, 5, 0.25)
    # x, y = np.meshgrid(x, y)
    # r = np.sqrt(x ** 2 + y ** 2)
    # z = np.sin(r)
    t = 10 ** -4.5

    u = 7.4e-08
    ch = 'chr1'
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(gmp, coords=[[10, 20]])
    bkc = BkgdCalculator(t, u, ano)

    tf = root_dir + '/precalc/littleB50Debug_2/blk.intg.lookup.table.txt'
    with open(tf, 'r') as f:
        f.readline()
        ar = []
        for line in f:
            line = line.split()
            x = float(line[0].split(',')[0].split('=')[1])

            al = []
            for li in line[1:]:
                li = li.split(',')
                # block length to test
                b_len = 1000.0
                # y is the distance to the end of segment in M
                y = float(li[0].split('=')[1])
                # assume the segment is 2bp long and calculate M/Mb:
                m_mb = y / b_len
                # m_mb = 10**-8.7
                # b_len = max(1.0, int(y / m_mb))
                # correct z values by multiplying by M/Mb and dividing by udel
                z = float(li[1].split('=')[1]) * (-u / m_mb)
                # calculated z value using integral method
                z_calc_int = bkc.bkgd_intgr(0.0, x, x+y) / m_mb
                # calculate z value using sum method assuming 2 bases
                z_calc_sum = bkc.bkgd_sum(0.0, x + np.arange(b_len)*m_mb)
                al.append((x, y, z, z_calc_int, z_calc_sum))

            ar.append(al)

    ar = np.array(ar)
    # farray = root_dir + '/lsm_python/temp_data/lkup.tbl.npy'
    # np.save(farray, ar)
    # ar = np.load(farray)

    x = np.log10(ar[1:, :, 0])
    y = np.log10(ar[1:, :, 1])
    z = np.exp(ar[1:, :, 2])
    z_int = np.exp(ar[1:, :, 3])
    z_sum = np.exp(ar[1:, :, 4])

    fig = plt.figure(figsize=(13, 7))
    plt.subplots_adjust(top=1.00, left=0.00, right=1.00, bottom=0.00)
    ax = fig.gca(projection='3d')
    surf1 = ax.plot_surface(x, y, z - z_int, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)
    ax.set_zlabel('B_table - B_integral_exact', labelpad=15)
    ax.set_ylabel('length (log10 Morgans)', labelpad=10)
    ax.set_xlabel('distance (log10 Morgans)', labelpad=10)
    fig.colorbar(surf1, shrink=0.5, aspect=5)

    fig = plt.figure(figsize=(13, 7))
    plt.subplots_adjust(top=1.00, left=0.00, right=1.00, bottom=0.00)
    ax = fig.gca(projection='3d')
    surf2 = ax.plot_surface(x, y, z - z_sum, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)
    ax.set_zlabel('B_table - B_sum_exact', labelpad=15)
    ax.set_ylabel('length (log10 Morgans)', labelpad=10)
    ax.set_xlabel('distance (log10 Morgans)', labelpad=10)
    fig.colorbar(surf2, shrink=0.5, aspect=5)

    plt.show()


def segment_length_properties():
    # cst = ChromStruct('chr1')
    # # slens = []
    # rlens = []
    # for ch in cst.chroms[-1::-1]:
    #     cst.chrom = ch
    #     gf = cst.gmap_files
    #     # gf = gf.replace('10kb/', '10kb/unfilled_AA_Maps/')
    #     af = cst.bs_targets.values()[0]
    #     gmp = GeneticMap(ch, gf, gscale=0.01)
    #     # gmp.close_gaps()
    #     asg = AnnoSegments(ch, af, gmp)
    #     # slens.append(asg.lengths)
    #     rlens.append(asg.rlengths)
    #
    # rlens = np.concatenate(rlens)
    # slens = np.concatenate(slens)
    # assert np.all(slens > 0)

    rout = root_dir + '/lsm_python/temp_data/rlens.npz'
    sout = root_dir + '/lsm_python/temp_data/slens.npz'

    # np.savez_compressed(rout, rlens=rlens)
    # np.savez_compressed(sout, slens=slens)
    nb = 100

    r = np.load(rout)['rlens']
    s = np.load(sout)['slens']
    si = s.argsort()
    s, r = s[si], r[si]
    step = len(s) / nb
    # idx = np.arange(0, len(s), step)
    idx = np.searchsorted(s, 2.5 * np.logspace(0, 3, 25))
    idx[-1] = len(s)

    # length in M
    rst = np.log10(r[r > 0].min())
    ren = np.log10(r.max())
    rbins = np.concatenate(([0], np.logspace(rst, ren, nb, True)))

    means, rhist, rates = [], [], []
    for (i, j) in zip(idx[:-1], idx[1:]):
        mean = np.log10(np.mean(s[i:j]))
        # rate = r[i:j] / s[i:j]
        # rcnt = np.histogram(rate, rbins)[0] / float(len(s))
        rcnt = np.histogram(r[i:j], rbins)[0] / float(len(s))
        means.append(mean)
        # rates.append(np.log10(rate.mean()))
        rhist.append(rcnt)

    # plt.figure(figsize=(13, 7))
    # plt.scatter(means, rates)
    # plt.show()

    z = np.array(rhist).T
    x, y = np.meshgrid(np.array(means), np.log10(rbins[1:]))

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    fig = plt.figure(figsize=(13, 7))
    plt.subplots_adjust(top=1.00, left=0.00, right=1.00, bottom=0.00)
    ax = fig.gca(projection='3d')

    # Plot the surfaces
    surf1 = ax.plot_surface(x, y, z, cstride=1, rstride=1,
                            cmap=cm.coolwarm, linewidth=0,
                            antialiased=False)

    # Customize the z axis
    ax.set_zlabel('fraction', labelpad=15)
    ax.set_ylabel('genetic length (log10 Morgans/bp)', labelpad=10)
    ax.set_xlabel('length (log10 bp)', labelpad=10)

    # Add a color bar which maps values to colors
    fig.colorbar(surf1, shrink=0.5, aspect=5)

    plt.show()


def errors_binned_by_exact_b(bdir, aut):
    b_bins = 100
    e_bins = 1000

    t = 10 ** -4.5
    eps = 0.01
    run_eps = 0.01
    rmax = '0.0001'  # default value
    # rmax = '{:.5f}'.format(1e-5)
    # rmax = '{:.6f}'.format(1e-6)
    tk = 'pr95.cleanrun'
    anno = 'primate_cons95_Segments'

    # sf = 'edge'
    sf = 'rand'

    save_figs = False

    cst = ChromStruct('chr1', tkn=tk, bdir=bdir, bs_annos=(anno,))
    f_str = cst.bxct_file(anno, t, sf)

    b_arrs = [np.load(f_str.replace('chr1', c)) for c in aut]
    b_arrs = np.vstack(b_arrs)
    log_bxct_min = np.log10(b_arrs[:, 2].min())
    bins = np.logspace(np.floor(log_bxct_min), 0, b_bins)

    # errs = bmap_err(b_arrs[:, 2], b_arrs[:, 1])
    errs = calculate_error(b_arrs[:, 2], b_arrs[:, 1], fabs=False)
    hl, ol = hist_error(eps, e_bins, [errs])
    names = [r'$\mathrm{interpolation\ \epsilon=%g}$' '\n'
             r'$\mathrm{recomb\ max=%s}$' % (run_eps, rmax)]
    if sf == 'edge':
        tok = bdir + ' exact evals '
    else:
        tok = bdir + ' random interp '
    plot_error_dist(t, eps, hl, ol, names, tok)

    if save_figs:
        pdir = 'murphy_et_al/figures/precision_tests/err_dist'
        fimg = '{}/{}/B.err.dist.{}.{}.png'.format(root_dir, pdir, bdir, sf)
        plt.savefig(fimg, dpi=256)
        plt.close()

#    ridx = np.where(abs(errs) > 2 * eps)[0]
#    if ridx.any():
#        cnt, bins = np.histogram(b_arrs[ridx, 2], bins)
#        cnt = cnt.astype(float) / len(ridx)
#        plt.figure(figsize=(8, 4))
#        plt.subplots_adjust(bottom=0.17)
#        plt.step(bins[1:], cnt)
#        txt = 'distribution of exact B where error is outside tolerance range'
#        plt.title(txt)
#        plt.xlabel('exact B')
#        plt.ylabel('fraction')
#        plt.xscale('log')
#
#        if save_figs:
#            fimg = '{}/{}/B.hist.{}.{}.png'.format(root_dir, pdir, bdir, sf)
#            plt.savefig(fimg, dpi=256)
#            plt.close()
        
    else:
        plt.show()

    return None


def haldane_vs_linear_r(coeff):
    b_bins = 100
    # coeff = 10 ** -4.5
    tk = 'pr95.cleanrun'

    anno = 'primate_cons95_Segments'

    aut = human_autosomes[4:]
    # aut = ('chr7', 'chr8') + human_autosomes[9:]
    # aut = ['chr22', 'chr21']

    # bdir = 'test_mod_2'
    bdir = 'bscale_100'
    # bdir = 'bscale_200'

    sf = 'haldane'

    cst = ChromStruct('chr1', tkn=tk, bdir=bdir, bs_annos=(anno,))
    f_str = cst.bxct_file(anno, coeff, sf)

    b_arrs = [np.load(f_str.replace('chr1', c)) for c in aut]
    b_arrs = np.vstack(b_arrs)
    # log_bxct_min = np.log10(b_arrs[:, 2].min())
    # bins = np.logspace(np.floor(log_bxct_min), 0, b_bins)

    # little b haldane - linear (additive error)
    errs = calculate_error(b_arrs[:, 2], b_arrs[:, 1], rel=False, fabs=False)
    plt.figure(figsize=(15, 8))
    plt.hist(errs, bins=b_bins)
    plt.title('additive difference in little b with Haldane correction ' +
              r'$\mathrm{(t=10^{%.1f})}$' % np.log10(coeff), fontsize=20)
    plt.xlabel('little b haldane r - little b linear r', fontsize=20)
    plt.xticks(fontsize=20)
    plt.ylabel('difference counts', fontsize=20)
    plt.yticks(fontsize=20)

    pdir = '{}/murphy_et_al/figures/precision_tests/haldane_r'.format(root_dir)
    fimg = '{}/haldane_minus_linear_t{:.8f}.png'.format(pdir, coeff)
    plt.savefig(fimg, dpi=256)
    # plt.show()

    # pdir = 'murphy_et_al/figures/precision_tests/err_dist'
    # fimg = '{}/{}/B.err.dist.{}.{}.png'.format(root_dir, pdir, bdir, sf)
    # plt.savefig(fimg, dpi=256)
    # plt.close()

    # plt.show()


def interp_errors():
    t = 10 ** -4.5
    eps = 0.02
    ch = 'chr21'

    # aut = human_autosomes
    # aut = ('chr7', 'chr8') + human_autosomes[9:]

    # bdir = 'increase_offs10x'
    # bdir = 'increase_offs10x'
    # bdir = 'rm_norec'
    bdir = 'fix_table'
    # bdir = 'littleB50Debug_2'
    # bdir = 'littleB5016x_tbl_splt1'

    # an = 'low_rec_rate'
    an = 'primate_cons95_Segments'
    # an = 'primate_cons95_Segments_split'
    # an = 'rm_norec'

    tk = 'pr95.cleanrun'
    # sf = 'edge'
    sf = 'rand'

    # set path params
    cst = ChromStruct(ch, tkn=tk, bdir=bdir, bs_annos=(an,))
    f_interp = cst.bxct_file(an, t, sf)

    # load random interp sample result and exact map for chrom
    x_interp, b_interp, b_exact = np.load(f_interp).T

    # get x positions for errors outside of tolerance
    errs = calculate_error(b_exact, b_interp)
    x_errs = x_interp[errs > 2 * eps]
    x_pass = x_interp[errs <= 2 * eps]

    # find the annotation closest to error & pass sites
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(gmp, cst.bs_target(an))
    i_r = np.searchsorted(ano.start, x_errs)
    i_p = np.searchsorted(ano.start, x_pass)

    # get the min distance to annotations in genetic map units
    gx_errs = gmp.interp_gpos(x_errs)
    gr = np.minimum(gx_errs-ano.rend[i_r-1], ano.rstart[i_r]-gx_errs)
    gx_pass = gmp.interp_gpos(x_pass)
    gp = np.minimum(gx_pass-ano.rend[i_p-1], ano.rstart[i_p]-gx_pass)

    # get exact points binning the interp error
    f_xct = cst.bkgd_file(an, t).replace('.bkgd', '.exactB.txt')
    x_xct = np.loadtxt(f_xct)[:, 0]
    rjdx = np.searchsorted(x_xct, x_errs)
    xct_i, xct_j = x_xct[rjdx-1], x_xct[rjdx]
    assert np.all((x_errs >= xct_i) & (x_errs <= xct_j))

    # count the number of conserved sites between xct points
    # get genetic distance between points with errors between them
    # gxct_i, gxct_j = gmp.interp_gpos(np.array([xct_i, xct_j]))
    # err_gdist = gxct_j - gxct_i

    # set bounds for histogram bins (exclude 0 for log scale)
    gmin = min(gp[gp > 0].min(), gr[gr > 0].min())
    gmax = max(gp.max(), gr.max())
    lg_gmin = np.floor(np.log10(gmin))
    lg_gmax = np.ceil(np.log10(gmax))
    lg_bins = np.logspace(lg_gmin, lg_gmax, 100)

    gi_cnt = np.histogram(gp, lg_bins)[0]
    g_cnt = np.histogram(gr, lg_bins)[0]

    plt.figure(figsize=(13, 7))
    plt.step(lg_bins[1:], gi_cnt / float(gi_cnt.sum()), label='in range')
    plt.step(lg_bins[1:], g_cnt / float(g_cnt.sum()), label='in range')
    plt.xscale('log')
    plt.legend()
    plt.show()


def cons_in_xct(ch, bdir, anno):
    t = 10 ** -4.5
    eps = 0.02
    tk = 'pr95.cleanrun'
    sf = 'rand'

    # set path params
    cst = ChromStruct(ch, tkn=tk, bdir=bdir, bs_annos=(anno,))
    f_interp = cst.bxct_file(anno, t, sf)

    # load random interp sample result and exact map for chrom
    x_interp, b_interp, b_exact = np.load(f_interp).T

    # get x positions for errors outside of tolerance
    errs = calculate_error(b_exact, b_interp)
    x_errs = x_interp[errs > 2 * eps]
    x_pass = x_interp[errs <= 2 * eps]

    # find the annotation closest to error & pass sites
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(gmp, cst.bs_target(anno))

    # get exact points binning the interp error
    f_xct = cst.bkgd_file(anno, t).replace('.bkgd', '.exactB.txt')
    x_xct = np.loadtxt(f_xct)[:, 0]

    # TEMP PATCH FOR INCOMPLETE EXACT TABLE
    x_errs = x_errs[x_errs < x_xct.max()]
    x_pass = x_pass[x_pass < x_xct.max()]

    i_r = np.searchsorted(x_xct, x_errs)
    i_p = np.searchsorted(x_xct, x_pass)

    xct_ri, xct_rj = x_xct[i_r-1], x_xct[i_r]
    xct_pi, xct_pj = x_xct[i_p-1], x_xct[i_p]

    assert np.all((x_errs >= xct_ri) & (x_errs <= xct_rj))
    # assert np.all((x_pass >= xct_pi) & (x_pass <= xct_pj))

    # count up the conserved segments within each exact point pair
    r_ncons = np.array([np.sum((ano.start >= i) & (ano.end <= j)) for
                        (i, j) in izip(xct_ri, xct_rj)])

    p_ncons = np.array([np.sum((ano.start >= i) & (ano.end <= j)) for
                        (i, j) in izip(xct_pi, xct_pj)])

    return r_ncons, p_ncons


def hist_cons_in_xct(r_ncons, p_ncons):
    plt.figure(figsize=(13, 7))
    _, bins, _ = plt.hist(r_ncons, bins=1000, label='out range', histtype='step',
                          normed=1, lw=1.5, cumulative=1)
    plt.hist(p_ncons, bins=1000, label='in range', histtype='step', normed=1,
             lw=1.5, cumulative=1)
    plt.xlabel('number of conserved segments between exact points', fontsize=16)
    plt.xticks(fontsize=16)
    plt.ylabel('fraction of interpolations', fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend()
    plt.show()


def gdist_to_cons(ch, bdir, anno):
    t = 10 ** -4.5
    eps = 0.02
    tk = 'pr95.cleanrun'
    sf = 'rand'

    # set path params
    cst = ChromStruct(ch, tkn=tk, bdir=bdir, bs_annos=(anno,))
    f_interp = cst.bxct_file(anno, t, sf)

    # load random interp sample result and exact map for chrom
    x_interp, b_interp, b_exact = np.load(f_interp).T

    # get x positions for errors outside of tolerance
    errs = calculate_error(b_exact, b_interp)
    x_errs = x_interp[errs > 2 * eps]
    x_pass = x_interp[errs <= 2 * eps]

    # find the annotation closest to error & pass sites
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(gmp, cst.bs_target(anno))
    i_r = np.searchsorted(ano.start, x_errs)
    i_p = np.searchsorted(ano.start, x_pass)

    # get the min distance to annotations in genetic map units
    gx_errs = gmp.interp_gpos(x_errs)
    gx_pass = gmp.interp_gpos(x_pass)

    r_dist = np.minimum(gx_errs - ano.rend[i_r - 1], ano.rstart[i_r] - gx_errs)
    p_dist = np.minimum(gx_pass - ano.rend[i_p - 1], ano.rstart[i_p] - gx_pass)

    return r_dist, p_dist


def hist_gdist2cons(rdist, pdist):
    # set bounds for histogram bins (exclude 0 for log scale)
    gmin = min(pdist[pdist > 0].min(), rdist[rdist > 0].min())
    gmax = max(pdist.max(), rdist.max())
    lg_gmin = np.floor(np.log10(gmin))
    lg_gmax = np.ceil(np.log10(gmax))
    lg_bins = np.logspace(lg_gmin, lg_gmax, 100)

    p_cnt = np.histogram(pdist, lg_bins)[0]
    r_cnt = np.histogram(rdist, lg_bins)[0]

    plt.figure(figsize=(13, 7))
    plt.step(lg_bins[1:], p_cnt / float(p_cnt.sum()), label='in range')
    plt.step(lg_bins[1:], r_cnt / float(r_cnt.sum()), label='out range')
    plt.xscale('log')
    plt.xlabel('distance to nearest conserved block (M)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.ylabel('fraction of interpolations', fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend()
    plt.show()


def errors_binned_between_cons(ch, bdir, eps):
    t = 10 ** -4.5
    tk = 'pr95.cleanrun'
    sf = 'rand'
    anno = 'primate_cons95_Segments'

    # init basic struct
    cst = ChromStruct(ch, tkn=tk, bdir=bdir, bs_annos=(anno,))

    # load gmap and annotation classes
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(gmp, cst.bs_target(anno))

    # load interp file, calc relative error
    f_interp = cst.bxct_file(anno, t, sf)
    x_interp, b_interp, b_exact = np.load(f_interp).T
    errs = calculate_error(b_exact, b_interp)

    # convert to genetic map positions
    x_interp = gmp.interp_gpos(x_interp)

    # segregate x positions for errors outside and inside of tolerance
    x_errs = x_interp[errs > 2 * eps]
    x_pass = x_interp[errs <= 2 * eps]

    # get the genetic distance between cons blocks
    mid_len = ano.rstart[1:] - ano.rend[:-1]
    # get midpoint from chrom start to first block, last block to chrom end
    mid_0 = [ano.rstart[0] / 2.0]
    mid_n = ano.rend[-1] + (gmp.gpos.max() - ano.rend[-1]) / 2.0
    # get first n-1 midpoints (excludes midpoint from last block to chrom end)
    mid_pts = np.concatenate((mid_0, ano.rend[:-1] + mid_len / 2.0))
    # stack midpoints and rcoords and reshape to 1D
    b_0 = np.column_stack((mid_pts, ano.rcoords)).reshape(3*ano.num_segs)
    # add gmap start=0, last midpoint and gmap end
    bins = np.concatenate(([0.0], b_0, [mid_n, gmp.gpos.max()]))

    # bbp = np.sum(gmp.interp_pos(bins[1:-1:3]) - gmp.interp_pos(bins[:-2:3]))
    # abp = np.sum(gmp.interp_pos(bins[2::3]) - gmp.interp_pos(bins[1:-1:3]))
    #
    # print '{}: below={:.3e}, above={:.3e}'.format(ch, bbp, abp)
    # return

    # check cons at expected indices and bins have expected length
    assert np.all(bins[2:-3:3] == ano.rstart)
    assert np.all(bins[3:-2:3] == ano.rend)
    assert bins.size == ano.num_segs * 3 + 3

    # histogram count errors and pass sites
    rcnt = np.histogram(x_errs, bins)[0]
    pcnt = np.histogram(x_pass, bins)[0]

    # there shouldnt be any sites in conserved blocks (but there are some)
    # assert np.sum(rcnt[2:-2:3]) + np.sum(pcnt[2:-2:3]) == 0

    # return pass_mid_1, pass_mid_2, pass_con, err_mid_1, err_mid_2, err_con
    cnts = np.vstack((pcnt[:-1:3], pcnt[1::3], rcnt[:-1:3], rcnt[1::3]))
    cn_cnts = np.vstack((pcnt[2:-2:3], rcnt[2:-2:3]))

    return cnts, cn_cnts


def summarize_midpoint_bincounts(bd, aut):
    # eps = 0.005
    eps = 0.01

    sample = 25000
    save_figs = False

    tdir = root_dir + '/lsm_python/temp_data'
    f_cnt_sums = '{}/{}.n_{}.cnt_sums.npy'.format(tdir, bd, sample)

    cstr = ChromStruct('chr1')
    cnt_sums = np.zeros(4)
#    aut = cstr.chroms[12:]
    for ch in aut:
        cnts, cn_cnts = errors_binned_between_cons(ch, bd, eps)
        cnt_sums += np.sum(cnts, axis=1)
    np.save(f_cnt_sums, cnt_sums)
    cnt_sums = np.load(f_cnt_sums)

    fracs = np.zeros(4)
    fracs[:] = cnt_sums
    fracs[:2] /= np.sum(fracs[:2])
    fracs[2:] /= np.sum(fracs[2:])
    hmax = fracs.max()

    plt.figure(figsize=(13, 7))
    ttl = 'sum of random neutral points below and above cons midpoints'
    ttl += r' ($\mathrm{\epsilon=%.1g}$)' % eps
    plt.title(ttl, fontsize=20)

    # check for significant differences in counts and label plot
    for (i, n) in enumerate(map(np.sum, (cnt_sums[:2], cnt_sums[2:]))):
        if i:
            shft = 2
            x_i, h_i, c_i = [2, 3], fracs[2:], cnt_sums[2:]
            lb, co = r'$\mathrm{error\ >\ 2\epsilon}$', 'darkorange'
        else:
            shft = 0
            x_i, h_i, c_i = [0, 1], fracs[:2], cnt_sums[:2]
            lb, co = r'$\mathrm{error\ \leqslant\ 2\epsilon}$', 'purple'

        # plot bars and label actual counts inside each bar
        plt.bar(x_i, h_i, align='center', label=lb, color=co)
        fmt = dict(ha='center', fontsize=18, color='white', fontweight='bold')
        plt.text(shft, 0.25, 'n = {:.0f}'.format(c_i[0]), **fmt)
        plt.text(shft+1, 0.25, 'n = {:.0f}'.format(c_i[1]), **fmt)

        # test whether error count for midpoint- bins is significantly
        # higher than error count for midpoint+ bins, which could indicate
        # that the search algorithm is more effective when moving towards
        # conserved blocks than away from them. assume errors are equally
        # likely to fall on either side side of the midpoint, so that error
        # count (k) for m+ and m- bins can be approximated by a binomial
        # random variable X~Bin(n, p) where p=0.5 and n=total errors
        p = 0.5
        k = cnt_sums[shft]

        # calculate mean, var and sd for Binom(n, p)
        mu = p*n
        var = n*p*(1.0-p)
        sd = np.sqrt(var)

        # calculate P(X >= k) to test for significantly greater m- error
        from scipy.stats import binom
        bnm = binom(n=n, p=p)
        p_high = bnm.sf(k-1)  # p(X > k-1)
        # calculate P(X <= k) to test for significantly lower m- error
        p_low = bnm.cdf(k)

        # use a normal approximation of the binomial and calculate Z scores
        # corresponding to P(X >= k) and P(X <= k)
        from scipy.stats import norm
        nrm = norm(loc=0, scale=1)
        z = (k - mu) / sd
        z_high = nrm.sf(z)
        z_low = nrm.cdf(z)

        # check the minimum p values from binomial and normal
        p_bnm = min(p_high, p_low)
        p_nrm = min(z_high, z_low)

        p_val = min(p_bnm, p_nrm)
        p_str = 'p = {:.3g}'.format(p_val)

        # label significance tests
        txt = '*** ' + p_str if p_val < 0.05 else p_str
        x = shft+0.5
        fmt = dict(ha='center', fontsize=24)
        plt.text(y=hmax+0.06, x=x, s=txt, **fmt)
        plt.axhline(y=hmax+0.03, xmin=x/4, xmax=(1+x)/4, lw=2.0, color='k')

    plt.axvline(x=1.5, lw=2.0, ls='--', color='dimgray')
    plt.xticks(range(4), ['below', 'above'] * 2, fontsize=20)
    plt.xlim(-0.5, 3.5)
    plt.ylabel('fraction of interpolations', fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(0, hmax+0.35)
    plt.legend(prop=dict(size=24), loc=9, ncol=2, columnspacing=8)

    if save_figs:
        pdir = 'murphy_et_al/figures/precision_tests/err_dist'
        fimg = '{}/{}/midpoint.err.{}.n_{}.png'.format(root_dir, pdir, bd, 
                sample)
        plt.savefig(fimg, dpi=256)
    else:
        print 'name: {}'.format(bd)
        plt.show()

    return None

    # a = [cons_in_xct(c, bd, an) for c in cstr.chroms[-2:]]
    # a = [gdist_to_cons(c, bd, an) for c in cstr.chroms]
    # ap = np.concatenate([x[1] for x in a])
    # ar = np.concatenate([x[0] for x in a])
    # np.save('{}/lsm_python/temp_data/rncons.npy'.format(cstr.root), ar)
    # np.save('{}/lsm_python/temp_data/pncons.npy'.format(cstr.root), ap)
    # ap = np.load('{}/lsm_python/temp_data/pncons.npy'.format(cstr.root))
    # ar = np.load('{}/lsm_python/temp_data/rncons.npy'.format(cstr.root))
    # hist_gdist2cons(ar, ap)
    # hist_cons_in_xct(ar, ap)


def main():
    s = 'std std_pts std_0.5prv_r std_pts_0.5prv_r newtab newtab_pts'
    aut = 'chr19 chr21 chr22'.split()
#    errors_binned_by_exact_b('std_haldane', aut)
    summarize_midpoint_bincounts('std_haldane', aut)
#    chroms = human_autosomes[6:10] + human_autosomes[11:]
#    for bd in s.split():
#        summarize_midpoint_bincounts(bd=bd, aut=chroms)
#        errors_binned_by_exact_b(bdir=bd, aut=chroms)
    # for t in [0.01, 0.001, 1e-4, 10**-4.5]:
    #     haldane_vs_linear_r(t)
    # plt.show()


if __name__ == '__main__':
    main()
    # main_1()
    # f = folders[0]
    # tk = tokens[0]
    # main_2(f, tk, 'chr1', 500)
    # main_3()
    # main_4()
    # main_5()
    # main_6()
    # segment_length_properties()
    # [idx2dbl(i) for i in xrange(250)]
