__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir
from figures.sub_rate_figs import neutcons_uratio_single
from skmisc import loess
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from figures.common_functions import format_panels, cst_from_fldr

style_dict = {'axes.facecolor': 'white', 'grid.color': 'k', 'grid.linestyle': '--', 'grid.linewidth': 0.2,
              'axes.edgecolor': 'k', 'axes.linewidth': 0.25, 'legend.fancybox': True, 'grid.alpha': 0.5,
              'lines.dashed_pattern': [10, 10], 'lines.scale_dashes': 1}
seaborn.set_theme(style_dict)
for (k, v) in style_dict.items():
    plt.rcParams[k] = v

# u0 for scaling parameter estimates
u0 = 1.4e-08

# make figures directory
final_dir = root_dir + '/result/final_files'
figdir = final_dir + '/mainfigs'
os.makedirs(figdir, exist_ok=True)


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def get_loess_line(fldr, span, return_points=False, load_con=False, loo=False):
    # load results in 2000 bins for LOESS plots
    fdir = final_dir + '/{}/'.format(fldr)
    if loo:
        sort_file = final_dir + '/{}/basic_sort_n{}_leaveOneOut.txt'.format(fldr, 2000)
    else:
        sort_file = final_dir + '/{}/basic_sort_n{}.txt'.format(fldr, 2000)
    div, pi, pred = np.loadtxt(sort_file).T
    # f_sort = fdir + 'sort_gc_cm_cn_il_n2000.txt'
    # div, pi, pred, gc, cn, cm, il = np.loadtxt(f_sort).T

    # load sorted conserved separately (done using 0.05cM radius of sites)
    if load_con:
        fcon = final_dir + '/{}/sorted.cons.n2000.txt'.format(fldr)
        cn = np.loadtxt(fcon)[:2000, 2]

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


#%% FIGURE 2A-B
def figure_2AB_final_format(fldr):
    """Figure 2A-B plot"""
    # style_dict = {'axes.facecolor': 'white', 'grid.color': 'k', 'grid.linestyle': '--', 'grid.linewidth': 0.2,
    #               'axes.edgecolor': 'k', 'axes.linewidth': 0.25, 'legend.fancybox': True, 'grid.alpha': 0.5,
    #               'lines.dashed_pattern': [10, 10], 'lines.scale_dashes': 1}
    # seaborn.set_theme(style_dict)
    # for (k, v) in style_dict.items():
    #     plt.rcParams[k] = v

    # fixed variables:
    color = 'darkorange'
    label = 'predicted (CADD)'
    clist = ['darkorange', 'dodgerblue']

    # set chr1 file name
    fchr1 = 'chr1.1.00e+06win_0.5slide_pred_and_obs.txt'
    f_name = final_dir + '/{}/{}'.format(fldr, fchr1)

    # load and normalize data
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    # remove some outliers
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan

    # create new figure
    plt.figure(figsize=(6.4, 2))
    plt.subplots_adjust(left=0.075, bottom=0.185, right=0.998, top=0.99,
                        wspace=0.32)

    # plot observed and predicted diversity in first 2/3 panel
    ax1 = plt.subplot(1, 3, (1, 2))

    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(xi, prd, color=color, lw=1.1, alpha=0.8, label=label)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    plt.yticks([0, 0.5, 1, 1.5], x=0.0125)
    # plt.ylim(0, 1.9)
    plt.ylim(0, 1.99)

    plt.xlabel('position on chromosome 1 (in Mb)', labelpad=2)
    xtck = range(25, 250, 25)
    xstr = [x if (x/25)%2 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    plt.xlim(0, xi[-1])
    # plt.tick_params(axis=u'both', which=u'both', length=0)

    legend = plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.2,
                        borderpad=0.2, framealpha=0.8, facecolor='white',
                        handlelength=1.2, labelspacing=0.25)
    legend.get_frame().set_linewidth(0.5)

    # plot R^2 in last 1/3 panel
    ax2 = plt.subplot(133)

    # load R^2 for just the current folder:
    rsq_folders = [fldr]
    # # load R^2 data for all conserved and exonic conserved
    # rsq_folders = [fldr, 'fish_cons94_new_jackknife_results',
    #                'cadd94_gmask_exonic', 'fish_cons94_new_exonic']
    rsq_shapes = ['o', '^', 'o', '^']
    rsq_labels = ['CADD', 'phastCons', r'$\mathrm{CADD_{e}}$',
                  r'$\mathrm{phastCons_{e}}$']
    rsq_colors = ['darkorange', 'deepskyblue', 'orangered', 'blue']
    rsq_list = []
    w = np.log10(np.loadtxt(final_dir + '/{}/rsq.log'.format(rsq_folders[-1]))[:16, 0])
    for fl in rsq_folders:
        if 'jackknife' in fl:
            fsq_file = final_dir + '/{}/{}.jkrsq.log'.format(fl, fl.replace('_jackknife_results', ''))
        else:
            fsq_file = final_dir + '/{}/rsq.log'.format(fl)
        rs = np.loadtxt(fsq_file)[:16, 1]
        rsq_list.append(rs)

    for i in range(len(rsq_folders)):
        rs = rsq_list[i]
        lab = rsq_labels[i]
        shape = rsq_shapes[i]
        col = rsq_colors[i]
        alpha = 0.6 if i in [1, 3] else 1.0
        plt.plot(w, rs, label=lab, color=col, marker=shape, lw=0, ms=5,
                 alpha=alpha)

    plt.xlabel('window size (log-scale)', labelpad=2)
    xtck = [4, 4.5, 5, 5.5, 6]
    xstr = [r'$10^{%.0f}$' % x if not x%1 else '' for x in xtck]
    plt.xticks(xtck, xstr, y=0.04)
    # plt.tick_params(axis=u'both', which=u'both', length=0)
    plt.ylim(0.05, 0.65)
    plt.ylabel(r'variance explained $(R^2)$', labelpad=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
    legend = plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
                        borderaxespad=0.2, ncol=1, frameon=1,
                        framealpha=0.8, facecolor='white', columnspacing=0.1,
                        labelspacing=0.25, prop=dict(size=8.5))
    legend.get_frame().set_linewidth(0.5)
    fig_letters = ['A', 'B']
    # fig_letters = ['a', 'b']

    plt.text(0.08, 0.915, fig_letters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.7475, 0.915, fig_letters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/fig2AB.{}.chr1.png'.format(fldr)
    plt.savefig(f_save, dpi=512)
    plt.close()


# set the results folder to look up results:
fldr = 'cadd94_gmask_v1.6_without_bstat'
figure_2AB_final_format(fldr)


#%% FIGURE 3
def figure_3(fldr, fcl, xlab, letter=None):
    # fixed variables:
    color = 'darkorange'
    label = 'predicted'

    # COLLATED PLOT
    f_col = final_dir + '/{}/{}.collate.5.00e-03width.npy'.format(fldr, fcl)
    bins, div, pi, pr, cnts = np.load(f_col).T
    obs = pi / cnts
    prd = pr / cnts
    y_max = obs.mean()
    obs /= y_max
    prd /= y_max
    bins += 2.5e-3

    plt.figure(figsize=(3.25, 1.95))
    plt.subplots_adjust(left=0.165, bottom=0.175, right=0.995, top=0.995)
    ax1 = plt.subplot(111)
    # format_panels(ax1)

    plt.plot(bins, obs, label='observed', color='darkslategray', lw=1.1)
    plt.plot(bins, prd, label=label, color=color, alpha=0.8, lw=1.1)
    plt.xlabel(xlab, labelpad=2)
    plt.xticks(y=0.05)
    # plt.xlim(-0.41, 0.41)
    plt.ylabel(r'diversity level ($\pi/\bar{\pi})$', labelpad=2)
    if fcl == 'YRI_syn_s1':
        plt.ylim(0.81, 1.13)
    if fcl == 'hc_derived_cons':
        plt.yticks([0.7, 1.0, 1.3], x=0.03)
    else:
        plt.yticks(x=0.03)
    legend = plt.legend(loc='lower right', ncol=1, frameon=1, borderaxespad=0.3,
                        borderpad=0.2, framealpha=0.75, facecolor='white',
                        handlelength=1.2, labelspacing=0.25)
    legend.get_frame().set_linewidth(0.5)
    if letter is not None:
        plt.text(0.145, 0.915, letter, transform=plt.gcf().transFigure)

    f_save = figdir + '/fig3.{}.{}.collated.png'.format(fldr, fcl)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat'
xlab = 'distance to nearest NS substitution (cM)'
figure_3(fldr, 'nonsyn', xlab)


#%% FIGURE 4 DATA
def figure_4_data():
    pct = [95, 94, 93]
    s_ape = [1 - neutcons_uratio_single('ape', pc=p)[0].mean() for p in pct]
    s_cad = [1 - neutcons_uratio_single('cadd', pc=p)[0].mean() for p in pct]
    s_fis = [1 - neutcons_uratio_single('fish', pc=p)[0].mean() for p in pct]

    # get inferred udel for ape and CADD
    u_fis = []
    u_cad = []
    u_ape = []
    for p in pct:
        if p == 94:
            cadfldr = 'cadd94_gmask_v1.6_without_bstat'
            fisfldr = 'fish_cons94_new'
            apefldr = 'ape_cons94_gmask_mnb_378'

        else:
            cadfldr = 'cadd{}_gmask_v1.6_without_bstat'.format(p)
            # cadfldr = 'cadd{}_gmask'.format(p)
            fisfldr = 'fish_cons{}_new'.format(p)
            apefldr = 'ape_cons{}_gmask'.format(p)

        u_cad.append(cst_from_fldr(cadfldr).stat.utot[0] * 1e8)
        u_fis.append(cst_from_fldr(fisfldr).stat.utot[0] * 1e8)
        u_ape.append(cst_from_fldr(apefldr).stat.utot[0] * 1e8)

    # return s_ape, s_cad, u_ape, u_cad
    return s_ape, s_cad, s_fis, u_ape, u_cad, u_fis


# FIGURE 4
def figure_4(s_fis, s_cad, u_fis, u_cad, changebackground):
    #  upper and lower bounds of u total
    ulo = 1.29
    uhi = 1.51
    pct = [95, 94, 93]

    # axis format variables for white background
    gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
    axis_lw = 0.25

    # plt.figure(figsize=(5, 3))
    fig = plt.figure(figsize=(3.25, 1.95))
    plt.subplots_adjust(wspace=0.05, left=0.2, right=1, top=0.995, bottom=0.21)
    # panel A: ape rates
    ax1 = plt.subplot(121)
    # version with white background
    if changebackground:
        ax1.set_facecolor('white')
        ax1.grid(**gridargs)
        ax1.axis('on')
        ax1.patch.set_edgecolor('black')
        ax1.patch.set_linewidth(axis_lw)

    x1 = np.arange(len(pct))
    utop = np.array([u / ulo for u in u_fis])
    ubot = np.array([u / uhi for u in u_fis])
    uheight = utop-ubot
    # print utop
    # print ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_fis, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('phastCons', labelpad=2)
    plt.ylabel('proportion of\ndeleterious mutations', labelpad=2)
    plt.yticks(x=0.06)
    plt.ylim(0, 1.15)

    # panel B: CADD rates
    ax2 = plt.subplot(122)
    # version with white background
    if changebackground:
        ax2.set_facecolor('white')
        ax2.grid(**gridargs)
        ax2.axis('on')
        ax2.patch.set_edgecolor('black')
        ax2.patch.set_linewidth(axis_lw)

    utop = np.array([u / ulo for u in u_cad])
    ubot = np.array([u / uhi for u in u_cad])
    uheight = utop-ubot
    lbl = 'background selection'
    plt.bar(x1, uheight, bottom=ubot, color='firebrick', alpha=0.75, label=lbl)
    lbl = 'evolutionary rates'
    plt.plot(x1, s_cad, 'o', color='k', ms=6, label=lbl)
    plt.plot([-0.5, len(pct)-0.5], [1, 1], 'k--', alpha=0.6)
    xtck = ['{}%'.format(int(100 - p)) for p in pct]
    plt.xticks(x1, xtck, rotation=45, y=0.065)
    plt.xlim(-0.45, len(pct)-0.55)
    plt.xlabel('CADD', labelpad=2)
    plt.yticks(color='none', x=0.05)
    plt.ylim(0, 1.15)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=6.5))
    # leg._legend_box.align = 'left'
    # plt.setp(leg.get_title(), fontsize='xx-small')
    # customize legend
    handles, labels = ax1.get_legend_handles_labels()
    leg = fig.legend(reversed(handles), reversed(labels), loc='lower center',
                     title='estimates from:', frameon=1, framealpha=1,
                     facecolor='white', prop=dict(size=8.5), ncol=2,
                     borderpad=0.2, labelspacing=0.2, handlelength=1.1,
                     handletextpad=0.4, bbox_to_anchor=(0.6, 0.195),
                     columnspacing=0.8)
    # leg = plt.legend(title='estimates from:', loc='lower left', frameon=1,
    #                  framealpha=0.75, facecolor='white', prop=dict(size=9),
    #                  ncol=2)
    leg._legend_box.align = 'left'
    leg.set_title('estimates from:', prop={'size':8.5})
    fig_letters = ['A', 'B']

    plt.text(0.21, 0.92, fig_letters[0], transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.62, 0.92, fig_letters[1], transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = figdir + '/fig4.urates.png'

    plt.savefig(f_save, dpi=512)
    plt.close()


# load the data for fig. 4
s_ape, s_cad, s_fis, u_ape, u_cad, u_fis = figure_4_data()
# create the figure
figure_4(s_fis, s_cad, u_fis, u_cad, changebackground=True)


#%% FIG. 5: OBSERVED VS. PREDICTED FINAL FORMATING
def figure_5_final_format(fldr, span=0.1, loo=False):
    """create loess smoothed plots for conservation and recombination rates"""
    # get 100 points data
    plt.rcParams['hatch.linewidth'] = 0.75
    if loo:
        sort_file = final_dir + '/{}/basic_sort_n{}_leaveOneOut.txt'.format(fldr, 100)
    else:
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

    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    # format_panels(ax1)

    axmin, axmax = 0.55, 1.19
    tx1, tx2 = 0.105, 0.61
    ty1, ty2 = 0.94, 0.48
    xtick = np.arange(0.5, 1.2, 0.1)

    # plot y=x line
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='k', ls='--', dashes=(5, 3))
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    # plot 5 bins affected by thresholding as hollow
    plt.scatter(pred[:5], pi[:5], marker='o', s=25, edgecolors='darkorange', linestyle='-',
                hatch='//'*6, alpha=0.4, color='none', linewidths=0.75)
    # plt.plot(pred[:5], pi[:5], marker='o', ms=5, color='gray', lw=0,
    #          alpha=0.5, markerfacecolor='None')
    # plot remaining 95 bins as solid
    plt.plot(pred[5:], pi[5:], marker='o', ms=5, color='darkorange', lw=0,
             alpha=0.5)
    # make a big dot where 5 bins are averaged
    plt.scatter(np.mean(pred[:5]), np.mean(pi[:5]), s=30, alpha=0.9, marker='o', color='darkred', zorder=3)
    # print(pi[0:10])
    print(pred[:5])
    print(pred[-5:])
    # plot LOESS line
    plt.plot(prlo, pilo, lw=2, color='orangered')

    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xticks(xtick, y=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0\ (B)$', labelpad=3)
    plt.ylim(0.45, axmax)
    plt.xlim(axmin, 1.02)
    # plt.xlim(0.45, 1.45)
    # solve y=x rotation
    adlen = axmax - 0.52
    oplen = adlen * (1.02 - axmin) / (axmax - 0.52)
    rot = np.arctan((oplen / adlen)) * (180.0 / np.pi)
    plt.text(0.75, 0.81, r'$y=x$', rotation=rot, ha='center', va='center',
             color='k')
    # plt.text(0.75, 0.81, r'$y=x$', rotation=38, ha='center', va='center',
    #          color='darkslategray', alpha=0.65)
    # plt.text(tx1, ty1, 'A', transform=plt.gcf().transFigure)
    if loo:
        f_save = figdir + '/fig5.{}.leave-one-out.png'.format(fldr)
    else:
        f_save = figdir + '/fig5.{}.png'.format(fldr)

    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_v1.6_without_bstat'
figure_5_final_format(fldr)


#%%