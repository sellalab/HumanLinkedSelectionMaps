__author__ = 'davidmurphy'


import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir, np
from figures.common_functions import hist_error, calculate_error, final_dir, \
    format_panels


#%% OLD VERSIONS OF PLOTTING FUNCTIONS
def plot_error_dist(t, eps, hlist, olist, names, token):
    tol = 2.0 * eps  # tolerance range = 2 x epsilon
    map_cols = ['fuchsia', 'darkorange', 'darkturquoise']
    bar_cols = np.array([
        ['fuchsia', 'blueviolet', 'darkviolet', 'purple', 'indigo'],
        ['darkorange', 'tomato', 'crimson', 'firebrick', 'darkred']])
    sdict = dict(linewidth=0.8, alpha=0.75)

    # plot error distribution
    # fig = plt.figure(figsize=(6.4, 4.8))
    fig = plt.figure(figsize=(3.25, 2.4))
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
    # plt.suptitle(title, fontsize=11)
    plt.suptitle(title, fontsize=11)

    # plot distribution of errors in leftmost 80% of plot area
    line_objs = []  # store line handles
    ax1 = fig.add_subplot(1, 5, (1, 4))
    format_panels(ax1)
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
    # plt.xticks(x_points, x_labels, fontsize=12)
    # plt.xlabel('relative error', fontsize=12)
    # plt.yticks(fontsize=11)
    # plt.ylabel('fraction', fontsize=12)
    # ax1.legend(line_objs, names, loc=2, fontsize=10)
    plt.xticks(x_points, x_labels)
    plt.xlabel('relative error')
    plt.ylabel('fraction')
    ax1.legend(line_objs, names, loc=2)

    # plot bar plots of outlier values stratified into epsilon * 2^x bins
    ax2 = fig.add_subplot(155)
    format_panels(ax2)
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
    # plt.xticks([0.5, ], ['{:.2g}%'.format(x) for x in fracs], fontsize=11)
    # plt.yticks(fontsize=11)
    # plt.xlim(0, 1)
    # plt.ylabel(r'$\mathrm{err>2\epsilon\ (\%)}$', fontsize=12)
    plt.xticks([0.5, ], ['{:.2g}%'.format(x) for x in fracs])
    plt.xlim(0, 1)
    plt.ylabel(r'$\mathrm{err>2\epsilon\ (\%)}$')
    plt.grid('off')
    # plt.show()


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

    save_figs = True

    cst = ChromStruct('chr1', tkn=tk, bdir=bdir, bs_annos=(anno,))
    f_str = cst.bxct_file(anno, t, sf)

    b_arrs = [np.load(f_str.replace('chr1', c)) for c in aut]
    b_arrs = np.vstack(b_arrs)
    log_bxct_min = np.log10(b_arrs[:, 2].min())
    bins = np.logspace(np.floor(log_bxct_min), 0, b_bins)

    # errs = bmap_err(b_arrs[:, 2], b_arrs[:, 1])
    errs = calculate_error(b_arrs[:, 2], b_arrs[:, 1], fabs=False)
    emask = (abs(errs)<=eps)
    errs = errs[emask]
    hl, ol = hist_error(eps, e_bins, [errs])
    names = [r'$\mathrm{interpolation\ \epsilon=%g}$' '\n'
             r'$\mathrm{recomb\ max=%s}$' % (run_eps, rmax)]
    if sf == 'edge':
        tok = bdir + ' exact evals '
    else:
        tok = bdir + ' random interp '
    plot_error_dist(t, eps, hl, ol, names, tok)

    if save_figs:
        fimg = final_dir + '/sfigs/B.err.dist.{}.{}.png'.format(bdir, sf)
        plt.savefig(fimg, dpi=512)
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


#%% BINNED ERRORS
def errors_binned_by_exact_b_2(bdir, aut):
    e_bins = 1000

    t = 10 ** -4.5
    eps = 0.01
    tk = 'pr95.cleanrun'
    anno = 'primate_cons95_Segments'

    if bdir == 'std_split_pts':
        sf = 'rand_merged'
    else:
        sf = 'rand'

    cst = ChromStruct('chr1', tkn=tk, bdir=bdir, bs_annos=(anno,))
    f_str = cst.bxct_file(anno, t, sf)

    b_arrs = [np.load(f_str.replace('chr1', c)) for c in aut]
    b_arrs = np.vstack(b_arrs)

    errs = calculate_error(b_arrs[:, 2], b_arrs[:, 1], fabs=False)
    if bdir == 'std_haldane_expm1_pts':
        emask = (abs(errs)<=eps)
        errs = errs[emask]
    hl, ol = hist_error(eps, e_bins, [errs])

    return hl, ol


def errors_binned_by_exact_b_3(bdir, aut):
    e_bins = 1000

    t = 10 ** -4.5
    eps = 0.01
    tk = 'pr95.cleanrun'
    anno = 'primate_cons95_Segments'

    if bdir == 'std_split_pts':
        sf = 'rand_merged'
        cst = ChromStruct('chr1', tkn=tk, bdir=bdir, bs_annos=(anno,))
        f_str = cst.bxct_file(anno, t, sf)

        b_arrs = [np.load(f_str.replace('chr1', c)) for c in aut]
        b_arrs = np.vstack(b_arrs)

        errs = calculate_error(b_arrs[:, 2], b_arrs[:, 1], fabs=False)
        if bdir == 'std_haldane_expm1_pts':
            emask = (abs(errs) <= eps)
            errs = errs[emask]
        hl, ol = hist_error(eps, e_bins, [errs])
        return hl, ol

    if bdir == 'exactPyB_vs_exactMcvB':
        f_str = root_dir + '/precalc/{}/chr1.PyVsC.0.000032.npy'.format(bdir)
        b_arrs = [np.load(f_str.replace('chr1', c)) for c in aut]
        b_arrs = np.vstack(b_arrs)
        errs = calculate_error(b_arrs[:, 2], b_arrs[:, 1], fabs=False)
        hl, ol = hist_error(eps, e_bins, [errs])
        return hl, ol
#
# # aut = 'chr19 chr21 chr22'.split()
# # bdir = 'std_haldane_expm1_pts'
# aut = ['chr22']
# bdir = 'increase_eps'
# hl, ol = errors_binned_by_exact_b_2(bdir, aut)
#

#%% PLOT ERROR DISTRIBUTION
def plot_error_dist_2(hlist, olist, token, title, letter):
    eps = 0.01
    tol = 2.0 * eps  # tolerance range = 2 x epsilon
    map_cols = ['fuchsia', 'darkorange', 'darkturquoise']
    bar_cols = np.array([
        ['fuchsia', 'blueviolet', 'darkviolet', 'purple', 'indigo'],
        ['darkorange', 'tomato', 'crimson', 'firebrick', 'darkred']])
    sdict = dict(linewidth=0.8, alpha=0.75)
    names = [r'$t=10^{-4.5}$']
    # plot error distribution
    # fig = plt.figure(figsize=(6.4, 4.8))
    fig = plt.figure(figsize=(3.25, 3))
    fig.subplots_adjust(left=0.2, bottom=0.1, right=0.87, top=0.93,
                        wspace=0.15)

    # plot distribution of errors in leftmost 80% of plot area
    line_objs = []  # store line handles
    ax1 = fig.add_subplot(1, 5, (1, 4))
    plt.title(title, y=0.99)
    format_panels(ax1)
    ymax = 0
    for (i, h) in enumerate(hlist):
        freq = h[0].astype(float) / h[0].sum()
        ymax = max(ymax, freq.max())
        bins = h[1][1:]
        line, = ax1.step(bins, freq, color=map_cols[i], **sdict)
        line_objs.append(line)
    # ax1.axvline(-tol, color='k', ls='--', alpha=0.5, lw=0.6)
    # xtol = ax1.axvline(tol, color='k', ls='--', alpha=0.5, lw=0.6)
    # names.append('jump size')
    # line_objs.append(xtol)

    ax1.axvline(-eps, color='k', ls=':', alpha=0.9, lw=0.8)
    xtol = ax1.axvline(eps, color='k', ls=':', alpha=0.9, lw=0.8)
    # names.append('precision bound')
    names.append(r'preset precision ($\epsilon$)')

    line_objs.append(xtol)
    # plot x labels in units of epsilon
    x_points = [x*tol for x in -2, -1, 0, 1, 2]
    # x_labels = [r'$\mathrm{%d\epsilon}$' % x if x != 0 else r'$\mathrm{0}$'
    #             for x in -4, -2, 0, 2, 4]
    x_labels = [-4, -2, 0, 2, 4]
    plt.xticks(x_points, x_labels, y=0.03)
    plt.xlabel(r'relative error (units of $\epsilon$)', labelpad=2)
    plt.yticks(x=0.04)
    plt.ylabel('fraction of test points', labelpad=2)
    plt.ylim(0, 1.25*ymax)
    ax1.legend(line_objs, names, loc=2, frameon=True, framealpha=0.75,
               facecolor='white', handlelength=1.2)

    # plot bar plots of outlier values stratified into epsilon * 2^x bins
    ax2 = fig.add_subplot(155)
    format_panels(ax2)
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
                txt = '>{:.0f}'.format(2**j) + r'$\mathrm{\epsilon}$'
                xi, yi = xs[i]-0.15, bots[i] + 0.5 * p
                plt.text(x=xi, y=yi, s=txt, fontsize=8, color='k')
                bots[i] += p

    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    plt.yticks(x=0.8)
    plt.xticks([0.5, ], ['{:.2g}%'.format(x) for x in fracs], y=0.03)
    plt.xlim(0, 1)
    plt.ylabel(r'$\mathrm{errors>\epsilon\ (\%)}$')
    plt.grid('off')
    # plt.show()
    plt.text(0.01, 0.945, letter, transform=plt.gcf().transFigure,
             fontweight='bold')
    # fimg = final_dir + '/sfigs/updatedApril2021.B_error_{}.png'.format(token)
    fimg = final_dir + '/sfigs/updatedJune2021.B_error_{}.png'.format(token)

    plt.savefig(fimg, dpi=512)
    plt.close()


# # # # aut = 'chr19 chr21 chr22'.split()
aut = ['chr22']
titles = ['McVicker et al. algorithm', 'Modified algorithm']
# # # bdirs = ['increase_eps', 'std_haldane_expm1_pts']
bdir = 'std_haldane_expm1_pts'
hl, ol = errors_binned_by_exact_b_2(bdir, aut)
plot_error_dist_2(hl, ol, bdir, titles[1], 'b')

# aut = 'chr19 chr21 chr22'.split()
# aut = ['chr{}'.format(c) for c in xrange(4, 23)]
aut = ['chr4']
bd = 'exactPyB_vs_exactMcvB'
hl, ol = errors_binned_by_exact_b_3(bd, aut)
plot_error_dist_2(hl, ol, bd, titles[0], 'a')
#%%