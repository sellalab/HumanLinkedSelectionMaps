__author__ = 'davidmurphy'

import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, cst_from_fldr
from data_processing.functions import rsquared_function
from scipy.stats import linregress, pearsonr
from skmisc import loess
from figures.common_functions import format_panels


final_dir = root_dir + '/result/final_files'
u0 = 1.4e-08
col_dict = {'AFR': 'darkorange', 'EUR': 'steelblue', 'AMR': 'fuchsia',
        'EAS': 'purple', 'SAS': 'darkturquoise'}


#%% FUNCTIONS USED IN ANALYSES
def get_pop_list():
    f_pops = root_dir + '/data/snps/allpops.txt'
    with open(f_pops, 'r') as f:
        pops = [line.strip('\n') for line in f.readlines()]

    return pops


def get_pop_dict():
    pfile = root_dir + '/data/snps/grouped_pops.txt'
    # create list of pops ordered by group and pop:group dict
    pdict = {}
    # pops, grps = [], []
    with open(pfile, 'r') as f:
        for line in f:
            pp, gp = line[:-1].split()
            pdict[pp] = gp

    return pdict


def residual_term(yi, fi):
    # cov(e,f)/var(f)
    # ei = yi-fi
    # return np.cov(ei, fi) / np.var(fi)
    ss_tot = np.sum(np.square(yi - yi.mean()))  # total sum of squares -- proportional to variance in xobs data
    return np.sum(2 * yi * fi - 2 * fi ** 2 + 2 * fi * yi.mean() - 2 * yi * yi.mean()) / ss_tot


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def outlier_mask(prd, lbound, ubound):
    """mask to remove data outside of some prediction bounds"""
    mask = (prd >= lbound) & (prd <= ubound)
    return mask


def fraction_removed(prd, lbound, ubound):
    """function that calculates the fraction of False sites in a mask"""
    mask = outlier_mask(prd, lbound, ubound)
    return 1.0 - (1.0 * np.sum(mask) / len(prd))


def get_pop_array(win, pop, anno, arr_type='YRI'):
    """get obs/pred arrays from self-map or YRI-map for a single pop"""
    # set file name depending on array type being used
    if arr_type == 'YRI':
        f_str = '/arrays/yri_pred_obs_{:.3e}.npz'.format(win)
    else:
        f_str = '/arrays/pred_obs_{:.3e}.npz'.format(win)

    if pop == 'YRI':
        if 'cadd' in anno:
            d_str = '/{}'.format(anno)
            fa = final_dir + d_str + f_str
        else:
            d_str = '/{}'.format(anno)
            fa = final_dir + d_str + f_str
    else:
        d_str = '/{}_{}'.format(pop, anno.replace('.', '_'))
        fa = final_dir + d_str + f_str

    prd, obs = np.load(fa)['prdobs'].T

    return obs, prd


def get_all_arrays(win, anno, arr_type='YRI'):
    """get obs/pred arrays from self-map or YRI-map for each pop"""
    obs_list, prd_list = [], []
    for pop in get_pop_list():
        obs, prd = get_pop_array(win, pop, anno, arr_type)
        obs_list.append(obs)
        prd_list.append(prd)

    return obs_list, prd_list


def get_rsq_from_both_maps(win, pop, anno):
    obs1, pred1 = get_pop_array(win, pop, anno, 'YRI')
    rsq1 = rsquared_function(obs1, pred1)
    obs2, pred2 = get_pop_array(win, pop, anno, 'self')
    rsq2 = rsquared_function(obs2, pred2)

    return rsq1, rsq2


def get_rsq_across_windows(anno, wins):
    win_dict = {}
    for win in wins:
        rsq_vals = []
        for pop in get_pop_list():
            rsq_vals.append(get_rsq_from_both_maps(win, pop, anno))
        win_dict[win] = np.array(rsq_vals)

    return win_dict


def get_self_beta(slf_data, pop_idx, trim_idx=None):
    # get beta for self-population map
    obs, prd = slf_data[0][pop_idx], slf_data[1][pop_idx]
    # normalize both sets of values
    yi = obs / obs.mean()
    fi = prd / prd.mean()

    # get pred sorting index
    fi_sidx = np.argsort(fi)
    # sort data by pred
    fi = fi[fi_sidx]
    yi = yi[fi_sidx]
    ei = yi - fi

    # use trimmed data
    if trim_idx:
        lidx, hidx = trim_idx
        yi = yi[lidx:hidx]
        fi = fi[lidx:hidx]
        ei = ei[lidx:hidx]
    res_term = residual_term(yi, fi)

    # do linear regression on pred/obs
    m, yint, rval, pval, stderr = linregress(fi, ei)

    return m, res_term



#%% PARAMS FIGURE (USED IN PAPER S6)
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

    plt.title('c', fontweight='bold', loc='left', y=0.975)
    plt.title('(i)', fontweight='bold', loc='center', y=0.975)
    ymax = 0
    # plt.text(0.08, 0.83, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.08, 0.83, '(solid = phastCons 6%)', transform=plt.gcf().transFigure)
    # plt.text(0.08, 0.73, '(hatched = CADD 6%)', transform=plt.gcf().transFigure)
    # hatches = ['.', '/', 'x']

    # solid bar results
    istop = len(pmf_list) / 2
    leg1lines = []
    for (i, df) in enumerate(pmf_list[:istop]):
        assert len(df) == 1
        lbl = llist[i]
        col = clist[i]
        df = df[0]
        ymax = max(df.max(), ymax)
        # if i%2:
        li = plt.bar(xi + s, df, w, label=lbl, color=col, align='edge')[0]
        leg1lines.append(li)
        s += w

    leg1 = plt.legend(leg1lines, llist[:istop], loc='upper left', ncol=ncol,
                      title='phastCons 6%', borderaxespad=0.3,borderpad=0.2,
                      handlelength=1.2, labelspacing=0.25)

    # hatched bar results
    leg2lines = []
    for (i, df) in enumerate(pmf_list[istop:]):
        assert len(df) == 1
        lbl = llist[i]
        col = clist[i]
        df = df[0]
        ymax = max(df.max(), ymax)
        li = plt.bar(xi + s, df, w, color=col, label=lbl, align='edge',
                hatch=hatch_pattern)[0]
        leg2lines.append(li)
        s += w
    plt.legend(leg2lines, llist[istop:], loc='upper right', ncol=ncol,
               title='CADD 6%', borderaxespad=0.3,borderpad=0.2,
               handlelength=1.2, labelspacing=0.25)
    plt.gca().add_artist(leg1)

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
            # if i % 2:
            if i < len(flist) / 2:
                plt.bar(0+s, udel, w, color=clist[li])
            else:
                plt.bar(0+s, udel, w, color=clist[li], hatch=hatch_pattern)
                # plt.bar(0+s, udel, w, color=shaded_col, alpha=alpha_val)

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
        if i < len(flist) / 2:
            plt.bar(0 + s, p0, w, color=clist[i])
        else:
            plt.bar(0 + s, p0, w, color=clist[i], hatch=hatch_pattern)
            # plt.bar(0 + s, p0, w, color=shaded_col, alpha=alpha_val)
        s += w
    plt.xticks([])

    f_save = final_dir + '/sfigs/updateJune2021.different.pop.params.png'
    plt.savefig(f_save, dpi=512)


pdict = get_pop_dict()
anno1 = 'fish_cons94_new'
anno2 = 'cadd94_gmask_v1_6_without_bstat'

select_pops = ['CEU', 'GIH', 'JPT', 'MXL']
# flist = ['{}_mnb_378'.format(anno1), '{}_mnb_378'.format(anno2)]
# llist = ['YRI', 'YRI']
# for pop in select_pops:
#     flist += ['{}_{}'.format(pop, anno1), '{}_{}'.format(pop, anno2)]
#     llist += [pop, pop]
# clist = [col_dict[pdict[pop]] for pop in llist]
flist = ['{}'.format(anno1)]
llist = ['YRI']
for pop in select_pops:
    flist.append('{}_{}'.format(pop, anno1))
    llist.append(pop)
flist.append('cadd94_gmask_v1.6_without_bstat')
llist.append('YRI')
for pop in select_pops:
    flist.append('{}_{}'.format(pop, anno2))
    llist.append(pop)
clist = [col_dict[pdict[pop]] for pop in llist]
compare_params(flist, llist, clist)


#%% PEARSON CORR DATA (USED IN PAPER S6)
def get_corr_dicts(wins, pops, anno):
    corr_dict = {}
    for pop in pops:
        corrs = []
        for w in wins:
            yri_map = get_pop_array(w, pop, anno)[1]
            pop_map = get_pop_array(w, pop, anno, 'self')[1]
            pear = pearsonr(yri_map, pop_map)[0]
            corrs.append(pear)
        corr_dict[pop] = corrs

    return corr_dict


windows = [1.56e4, 3.12e4, 6.25e4, 1.25e5, 2.5e5, 5e5, 1e6]
pops = ['CEU', 'GIH', 'JPT', 'MXL']
corrdict1 = get_corr_dicts(windows, pops, 'cadd94_gmask_v1.6_without_bstat')
corrdict2 = get_corr_dicts(windows, pops, 'fish_cons94_new')


#%% PEARSON CORR FIGURE (USED IN PAPER S6)
def pearson_plots(peardict1, peardict2):
    """get the pearson correlations from each set of results and plot A & B"""
    windows = [1.56e4, 3.12e4, 6.25e4, 1.25e5, 2.5e5, 5e5, 1e6]
    use_pop = ['CEU', 'GIH', 'JPT', 'MXL']
    cols = ['steelblue', 'fuchsia', 'purple', 'darkturquoise']
    plt.figure(figsize=(6.5, 2.16))
    plt.rc('hatch', color='k', linewidth=0.25)
    plt.rc('axes.formatter', useoffset=False)
    plt.subplots_adjust(wspace=0.37, top=0.91, right=0.995, left=0.125,
                        bottom=0.25)
    ax1 = plt.subplot(121)
    format_panels(ax1)
    plt.title('phastCons 6%', y=0.97)
    xi = np.arange(len(windows))
    width = 0.8 / len(use_pop)
    offset = -0.4
    for (i, pop) in enumerate(use_pop):
        plt.bar(xi+offset, peardict1[pop], width, label=pop, color=cols[i])
        offset += width
    # plt.yticks([0.99992, 0.99994, 0.99996, 0.99998, 1], x=0.04)
    plt.ylabel('Pearson correlation')
    plt.ylim(0.999, 1.0001)
    xtck = ['10kb', '30kb', '60kb', '125kb', '250kb', '500kb', '1Mb']
    plt.xticks(xi, xtck, rotation=45, y=0.04)
    plt.xlabel('window size', labelpad=2)
    # plt.legend(loc='upper right', prop=dict(size=8.5), ncol=4, borderpad=0.2,
    #            labelspacing=0.2, handlelength=0.8, handletextpad=0.4,
    #            columnspacing=0.8, frameon=1, framealpha=0.75, facecolor='white')

    ax2 = plt.subplot(122)
    format_panels(ax2)
    plt.title('CADD 6%', y=0.97)
    offset = -0.4
    for (i, pop) in enumerate(use_pop):
        plt.bar(xi + offset, peardict2[pop], width, label=pop, color=cols[i],
                hatch='\\'*6)
        offset += width
    # plt.yticks([0.99992, 0.99994, 0.99996, 0.99998, 1], x=0.04)
    plt.ylabel('Pearson correlation')
    plt.ylim(0.999, 1.0001)
    xtck = ['10kb', '30kb', '60kb', '125kb', '250kb', '500kb', '1Mb']
    plt.xticks(xi, xtck, rotation=45, y=0.04)
    plt.xlabel('window size', labelpad=2)

    ytxt = 0.84
    plt.text(0.13, ytxt, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')
    plt.text(0.63, ytxt, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/updateJune2021.population.rsq.pearsons2.png'
    plt.savefig(f_save, dpi=512)


pearson_plots(corrdict2, corrdict1)


#%% RSQ VS YRI-RSQ PLOTS (USED IN PAPER S6)
def rsq_vs_yrirsq_plot(win_dict):
    windows = [1.1e4, 1.25e5, 1e6]
    plt.figure(figsize=(6.5, 2.16))
    plt.rc('axes.formatter', useoffset=False)
    plt.subplots_adjust(wspace=0.2, top=0.89, right=0.995, left=0.09,
                        bottom=0.16)
    pop_dict = get_pop_dict()
    colors = [col_dict[pop_dict[pop]] for pop in get_pop_list()]
    tcks = [[0.1, 0.12, 0.14, 0.16], [0.25, 0.3, 0.35, 0.4],
            [0.45, 0.5, 0.55, 0.6]]
    lims = [(0.095, 0.165), (0.26, 0.42), (0.44, 0.61)]
    ttls = ['10 kb windows', '100 kb windows', '1 Mb windows']
    seen_labs = []
    for i in xrange(3):
        ax = plt.subplot(1, 3, i+1)
        format_panels(ax)
        r1, r2 = win_dict[windows[i]].T
        plt.title(ttls[i], y=0.98)
        plt.plot([0, 1], [0, 1], lw=0.5, ls='--', color='darkslategray')
        for j in xrange(len(r1)):
            lab = pop_dict[get_pop_list()[j]]
            if lab not in seen_labs:
                seen_labs.append(lab)
                plt.plot(r1[j], r2[j], color=colors[j], alpha=0.75, label=lab,
                         marker='o', ms=6, lw=0)
            else:
                plt.plot(r1[j], r2[j], color=colors[j], alpha=0.75, marker='o',
                         ms=6, lw=0)

        # plt.scatter(r1, r2, color=colors, alpha=0.75)
        plt.xticks(tcks[i], y=0.04)
        plt.xlim(*lims[i])
        plt.xlabel(r'$R^2$ based on YRI map', labelpad=2)
        plt.yticks(tcks[i], x=0.04)
        plt.ylim(*lims[i])
        if i == 0:
            plt.ylabel(r'$R^2$ based on pop. map', labelpad=2)
            plt.legend(loc='upper left', prop=dict(size=8.5), ncol=1, frameon=1,
                       framealpha=0.8, facecolor='white', handletextpad=0,
                       borderaxespad=0.4)
    # plt.text(0.13, 0.92, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.64, 0.92, 'B', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/population.rsq.vs.rsq.png'
    plt.savefig(f_save, dpi=512)


win_dict = get_rsq_across_windows('cadd94_gmask_v1.6_without_bstat',
                                  [1.1e4, 1.25e5, 1e6])
rsq_vs_yrirsq_plot(win_dict)


#%% RSQ VS 1/VAR PLOTS (USED N PAPER S6)
def rsq_vs_var(win_dict, smap='YRI'):
    windows = [1.1e4, 1.25e5, 1e6]
    plt.figure(figsize=(6.5, 2.16))
    plt.rc('axes.formatter', useoffset=False)
    plt.subplots_adjust(wspace=0.2, top=0.89, right=0.995, left=0.09,
                        bottom=0.16)
    pop_dict = get_pop_dict()
    pop_list = get_pop_list()
    colors = [col_dict[pop_dict[pop]] for pop in pop_list]
    xtcks = [[2, 2.5, 3], [5, 6, 7, 8], [12, 14, 16, 18]]
    ytcks = [[0.10, 0.12, 0.14, 0.16], [0.25, 0.3, 0.35, 0.4], [0.45, 0.50, 0.55, 0.6, 0.65, 0.7]]
    ylims = [(0.095, 0.165), (0.26, 0.42), (0.44, 0.61)]
    xlims = [(1.6, 3.1), (4.7, 8.5), (11.75, 18.8)]
    xtxt = [0.095, 0.415, 0.735]
    ytxt = 0.82
    letters = ['a', 'b', 'c']
    ttls = ['10 kb windows', '100 kb windows', '1 Mb windows']
    # line_txt = r'$R^2=\frac{V(f)}{V(y)}(1+\beta)$'
    line_txt = r'$R^2=V(f)/V(y)$'
    seen_labs = []
    for i in xrange(3):
        ax = plt.subplot(1, 3, i+1)
        format_panels(ax)
        obs_list, prd_list = win_dict[windows[i]]
        plt.title(ttls[i], y=0.965, loc='center')
        # plt.text(xtxt[i], 0.81, letters[i], transform=plt.gcf().transFigure)
        # plt.plot([0, 1], [0, 1], lw=0.5, ls='--', color='darkslategray')
        for j in xrange(len(obs_list)):
            obs, prd = obs_list[j], prd_list[j]
            rsq = rsquared_function(obs, prd)
            var_obs = np.var(obs/obs.mean())
            var_prd = np.var(prd/prd.mean())
            # print var_prd
            # res = residual_term(obs/obs.mean(), prd/prd.mean())
            if j == 0:
                x1, x2 = xlims[i]
                y1, y2 = var_prd*x1, var_prd*x2
                ylims[i] = (y1, y2)
                plt.plot([x1, x2], [y1, y2], ls='--', lw=0.5,
                         color='darkslategray')
                if i == 0:
                    ymid = y1 + (y2-y1)/2.0
                    xmid = x1 + (x2-x1)/2.0
                    plt.text(1.02*xmid, 1.09*ymid, line_txt,
                             color='darkslategray', rotation=45, ha='center',
                             va='center', fontsize=10)
            lab = pop_dict[pop_list[j]]
            if lab not in seen_labs:
                seen_labs.append(lab)
                plt.plot(1.0/var_obs, rsq, color=colors[j], alpha=0.75, label=lab,
                         marker='o', ms=6, lw=0)
            else:
                plt.plot(1.0/var_obs, rsq, color=colors[j], alpha=0.75, marker='o',
                         ms=6, lw=0)
            # plt.plot(1.0/var_obs, rsq-res, marker='s', markerfacecolor='None',
            #  markeredgecolor=colors[j], markeredgewidth=0.9, ms=4)
            plt.text(xtxt[i], ytxt, letters[i],
                     transform=plt.gcf().transFigure,fontweight='bold')
        plt.xticks(xtcks[i], y=0.04)
        plt.xlim(*xlims[i])
        plt.xlabel(r'$1/V(y)$', labelpad=2)
        plt.yticks(ytcks[i], x=0.04)
        plt.ylim(*ylims[i])
        if i == 0:
            plt.ylabel(r'$R^2$', labelpad=2)
            plt.legend(loc='lower right', prop=dict(size=8.5), ncol=1, frameon=1,
                       framealpha=0.8, facecolor='white', handletextpad=0,
                       borderaxespad=0.4)
    # plt.text(0.13, 0.92, 'A', transform=plt.gcf().transFigure)
    # plt.text(0.64, 0.92, 'B', transform=plt.gcf().transFigure)

    f_save = final_dir + '/sfigs/updateJune2021.population.rsq.vs.var.{}.png'.format(smap)
    plt.savefig(f_save, dpi=512)
    plt.close()


win_dict_yri = {}
win_dict_slf = {}
for win in [1.1e4, 1.25e5, 1e6]:
    win_dict_yri[win] = get_all_arrays(win, 'cadd94_gmask_v1.6_without_bstat', 'YRI')
    win_dict_slf[win] = get_all_arrays(win, 'cadd94_gmask_v1.6_without_bstat', 'self')

rsq_vs_var(win_dict_yri)
# rsq_vs_var(win_dict_slf, 'self')


#%% SCATTER PLOTS (USED IN PAPER S6)

win = 1e6
yri_list = get_all_arrays(win, 'cadd94_gmask_v1.6_without_bstat', 'YRI')
slf_list = get_all_arrays(win, 'cadd94_gmask_v1.6_without_bstat', 'self')


#%% SCATTER PLOTS (USED IN PAPER S6)
def res_vs_pred(yri_data, slf_data, pop, letter):
    """scatter plot of predicted vs. residuals"""
    # get the index to the population being used
    pop_idx = get_pop_list().index(pop)
    # get the subset of values for the needed population from the YRI result
    obs, prd = yri_data[0][pop_idx], yri_data[1][pop_idx]
    # normalize both sets of values
    yi = obs / obs.mean()
    fi = prd / prd.mean()

    # print np.sum(np.square(yi-yi.mean()))

    # get pred sorting index
    fi_sidx = np.argsort(fi)
    # sort data by pred
    fi = fi[fi_sidx]
    yi = yi[fi_sidx]
    ei = yi-fi

    # get indices for lowest/highest 1.5% of data
    lidx = int(0.015*len(fi))
    hidx = int(0.985*len(fi))

    # do linear regression on pred/obs
    m, yint, rval, pval, stderr = linregress(fi, ei)

    # get the beta for the pop-specific maps
    slf_beta, slf_res = get_self_beta(slf_data, pop_idx)
    # print pop, m, slf_beta, residual_term(yi, fi), slf_res
    print pop, m, residual_term(yi, fi), yint

    # # do LOESS smoothing on pred/obs
    # wts = np.ones(shape=len(fi))
    # p_loess = predict_loess(fi, ei, wts, 0.1, fi)

    plt.figure(figsize=(2.16, 2.16))
    plt.subplots_adjust(top=0.995, right=0.995, left=0.18,
                        bottom=0.16)
    ax = plt.subplot(111)
    format_panels(ax)
    # plt.title(pop, y=0.87)
    # plot scatter of pred vs. residuals
    plt.plot(fi, ei, marker='o', ms=2.5,
             markerfacecolor='None', markeredgecolor='k',
             markeredgewidth=0.3, lw=0, alpha=0.75)

    # add shading to bottom/top 1.5% of windows
    outlier_col = 'deepskyblue'
    plt.plot(fi[:lidx], ei[:lidx], marker='o', ms=2.5,
             markerfacecolor='None', markeredgecolor=outlier_col,
             markeredgewidth=0.4, lw=0, alpha=1)
    plt.plot(fi[hidx:], ei[hidx:], marker='o', ms=2.5,
             markerfacecolor='None', markeredgecolor=outlier_col,
             markeredgewidth=0.4, lw=0, alpha=1)

    # plot regression line
    beta_lbl = r'$\beta_a=%.4f\ (%.4f)$' %(m, slf_beta)
    if pop == 'YRI':
        beta_lbl = r'$\beta_a=%.4f$' % m
    li = np.array([fi.min(), fi.max()])
    plt.plot(li, li*m+yint, color='red', ls='-', lw=0.75,
             label=beta_lbl)

    # repeat the analysis with points trimmed from extremes
    # do linear regression on pred/obs
    m, yint, rval, pval, stderr = linregress(fi[lidx:hidx], ei[lidx:hidx])

    # get the beta for the pop-specific maps
    slf_beta, slf_res = get_self_beta(slf_data, pop_idx, [lidx, hidx])
    # print pop, m, slf_beta, residual_term(yi[lidx:hidx], fi[lidx:hidx]), slf_res
    # plot regression line
    beta_lbl = r'$\beta_t=%.4f\ (%.4f)$' % (m, slf_beta)
    if pop == 'YRI':
        beta_lbl = r'$\beta_t=%.4f$' % m

    # li = np.array([fi[lidx], fi[hidx]])
    li = np.array([fi.min(), fi.max()])
    plt.plot(li, li * m + yint, color='magenta', ls=':', lw=0.75,
             label=beta_lbl)

    # # plot LOESS line
    # plt.plot(fi, p_loess, 'darkturquoise', lw=0.75, label='LOESS fit')

    # # plot moving average line
    # lbl += '_w{:.2f}_s{:.0f}'.format(wsize, 100*step)
    # slbl = 'sliding window\nwin={:.2f}; slide={:.0f}%'.format(wsize, 100.0*step)
    # plt.plot(mean_fi, mean_d_yf, 'darkturquoise', lw=0.75, marker='o',
    #          label=slbl, ms=3)

    plt.xticks(y=0.04)
    plt.xlabel(r'$f_i$ (1Mb windows)', labelpad=1)
    plt.yticks(x=0.04)
    plt.ylabel(r'$e_i$ (1Mb windows)', labelpad=2)
    plt.ylim(-1.6, 1.6)
    plt.xlim(0.51, 1.9)
    plt.legend(loc='lower center', prop=dict(size=8), handlelength=0.8,
               handletextpad=0.4, frameon=1, framealpha=0.75, facecolor='white')
    # plt.text(0.2, 0.92, letter, transform=plt.gcf().transFigure)
    # txt_str = '{}) {}'.format(letter, pop)
    # txt_str = '{}) {}'.format(letter, pop)
    plt.title(pop, y=0.92, x=0.61, transform=plt.gcf().transFigure)
    plt.text(0.2, 0.92, letter, transform=plt.gcf().transFigure,
             fontweight='bold')

    f_save = final_dir + '/sfigs/updateMarch2021.{}_res.vs.pred.png'.format(pop)
    plt.savefig(f_save, dpi=512)
    plt.close()


pops = ['YRI', 'CEU', 'GIH', 'JPT', 'MXL']
lets = 'd e f g h'.split()
for i in xrange(5):
    res_vs_pred(yri_list, slf_list, pops[i], lets[i])


#%% VAR(Y) VS. MEAN(PI) PLOT
def var_vs_mean(win_dict):
    windows = [1.1e4, 1.25e5, 1e6]
    plt.figure(figsize=(2.16, 2.16))
    plt.rc('axes.formatter', useoffset=False)
    plt.subplots_adjust(wspace=0.2, top=0.995, right=0.995, left=0.2,
                        bottom=0.16)
    pop_dict = get_pop_dict()
    pop_list = get_pop_list()
    colors = [col_dict[pop_dict[pop]] for pop in pop_list]
    seen_labs = []
    ax = plt.subplot(111)
    format_panels(ax)
    obs_list, prd_list = win_dict[1e6]
    for j in xrange(len(obs_list)):
        obs, prd = obs_list[j], prd_list[j]
        var_y = np.var(obs/obs.mean())
        mean_pi = np.mean(obs) * 1000
        lab = pop_dict[pop_list[j]]
        if lab not in seen_labs:
            seen_labs.append(lab)
            plt.plot(var_y, mean_pi, color=colors[j], alpha=0.75, label=lab,
                     marker='o', ms=6, lw=0)
        else:
            plt.plot(var_y, mean_pi, color=colors[j], alpha=0.75,
                     marker='o', ms=6, lw=0)

    plt.xticks([0.06, 0.07, 0.08], y=0.04)
    plt.xlabel('V(y)', labelpad=2)
    plt.yticks(x=0.04)
    plt.ylabel(r'$\bar{\pi}\ (\times10^3)$', labelpad=2)
    plt.legend(loc='upper right', prop=dict(size=8.5), ncol=1, frameon=1,
               framealpha=0.8, facecolor='white', handletextpad=0,
               borderaxespad=0.4)
    f_save = final_dir + '/sfigs/population.var.vs.meanpi.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


var_vs_mean(win_dict_yri)


#%%
def rsq_vs_mean(win_dict, smap='YRI'):
    plt.figure(figsize=(2.16, 2.16))
    plt.rc('axes.formatter', useoffset=False)
    plt.subplots_adjust(wspace=0.2, top=0.995, right=0.995, left=0.235,
                        bottom=0.16)
    pop_dict = get_pop_dict()
    pop_list = get_pop_list()
    colors = [col_dict[pop_dict[pop]] for pop in pop_list]
    seen_labs = []
    ax = plt.subplot(111)
    format_panels(ax)
    obs_list, prd_list = win_dict[1e6]
    r_list = []
    p_list = []
    for j in xrange(len(obs_list)):
        obs, prd = obs_list[j], prd_list[j]
        rsq = rsquared_function(obs, prd)
        mean_pi = np.mean(obs) * 1000
        lab = pop_dict[pop_list[j]]
        r_list.append(rsq)
        p_list.append(mean_pi)
        if lab not in seen_labs:
            seen_labs.append(lab)
            plt.plot(mean_pi, rsq, color=colors[j], alpha=0.75, label=lab,
                     marker='o', ms=6, lw=0)
        else:
            plt.plot(mean_pi, rsq, color=colors[j], alpha=0.75,
                     marker='o', ms=6, lw=0)

    # print linregress([0.001*p for p in p_list], r_list)
    m, yint, rval, pval, stderr = linregress(p_list, r_list)
    x1, x2 = min(p_list), max(p_list)
    y1, y2 = m*x1+yint, m*x2+yint
    plt.plot([x1, x2], [y1, y2], ls='--', color='darkslategray', lw=0.5)
    txt = r'$y=%.2fx+%.2f$' %(m, yint)
    plt.text(0.5, 0.32, txt, transform=plt.gcf().transFigure, va='bottom',
             ha='left', fontsize=10, rotation=49, color='darkslategray')
    plt.yticks(x=0.04)
    plt.ylabel(r'$\mathrm{R^2}$', labelpad=2)
    plt.ylim(0.44, 0.61)
    plt.xticks(y=0.04)
    plt.xlabel(r'$\bar{\pi}\ (\times10^3)$', labelpad=2)
    plt.legend(loc='upper left', prop=dict(size=8.5), ncol=1, frameon=1,
               framealpha=0.8, facecolor='white', handletextpad=0,
               borderaxespad=0.4)
    f_save = final_dir + '/sfigs/population.rsq.vs.meanpi.{}.png'.format(smap)
    plt.savefig(f_save, dpi=512)
    plt.close()


rsq_vs_mean(win_dict_yri)
# rsq_vs_mean(win_dict_slf, 'self')


#%%
# wstr = '1Mb'
# rsqlist, varlist, rsdlist = [], [], []
# nprd_list = []
# nobs_list = []
# pct_mask = []
# use_mask = 0
# lower_bound = 0.65
# upper_bound = 1.5
#
# for pop, prd, obs in zip(pops, prd_list, obs_list):
#     # normalize predicted and observed
#     nprd = prd / prd.mean()
#     nobs = obs / obs.mean()
#
#     # get mask for lower and upper bounds is specified
#     if use_mask:
#         # msk = abs(nprd-nobs) <= 0.4
#         sidx_nprd = np.argsort(nprd)
#         sort_nprd = nprd[sidx_nprd]
#         lowval = sort_nprd[int(0.025*len(nprd))]
#         highval = sort_nprd[int(0.975*len(nprd))]
#         msk = outlier_mask(nprd, lowval, highval)
#         # msk = outlier_mask(nprd, lower_bound, upper_bound)
#         nprd = nprd[msk]
#         nobs = nobs[msk]
#         nprd /= nprd.mean()
#         obs /= nobs.mean()
#         pmsk = 1.0 - (1.0 * np.sum(msk) / len(msk))
#     else:
#         pmsk = 0
#
#     # append masked (or unmasked) normalized values
#     pct_mask.append(pmsk)
#     nprd_list.append(nprd)
#     nobs_list.append(nobs)
#
#     # calculate variance, R^2 and residual term of normalized diversity
#     var = np.var(nobs)
#     rsq = rsquared_function(nobs, nprd)
#     rsd = residual_term(nobs, nprd)
#
#     rsqlist.append(rsq)
#     rsdlist.append(rsd)
#     varlist.append(var)
#     # print pop, rsq, (1.0/(cov**2)) * np.var(nprd), rsd
#
# rsqlist = np.array(rsqlist)
# rsdlist = np.array(rsdlist)
# varlist = np.array(varlist)
# nprd_list = np.array(nprd_list)
# nobs_list = np.array(nobs_list)
# pvar = np.var(nprd_list[0])
# assert np.allclose(pvar, np.var(nprd_list, axis=1))
#
#
# # RESIDUALS PLOT
# def cross_term_correction_plot(vi, rsq, rsd, pred_var, sname, remove=0):
#     plt.figure(figsize=(2.5, 2.5))
#     plt.subplots_adjust(top=0.98, right=1, bottom=0.2, left=0.24)
#     win_lab = '({} windows)'.format(wstr)
#     if remove > 0:
#         rlab = r'$R^2$ ({:.2f}% removed)'.format(remove)
#     else:
#         rlab = r'$R^2$'
#     plt.plot(1/vi, rsq, color='dodgerblue', label=rlab, marker='o',
#              markersize=5, lw=0)
#     plt.plot(1/vi, rsq-rsd, color='darkorange', marker='o',
#              markersize=5, lw=0, label=r'$R^2-cross\ term$')
#     # xi = np.array([11.5, 25.0])
#     xmin, xmax = 4.2, 10
#     ymin, ymax = 0.18, 0.42
#     # xmin, xmax = min(1/vi), max(1/vi)
#     xi = np.array([xmin, xmax])
#     plt.plot(xi, xi*pred_var, ls='--', color='gray', label='y=Var(B)x')
#     plt.ylabel(r'$R^2$ '+ win_lab)
#     plt.xlabel(r'$1/Var(\pi)$ ' +  win_lab)
#     plt.xlim(xmin, xmax)
#     plt.ylim(ymin, ymax)
#     plt.legend(frameon=1, framealpha=0.75, facecolor='white')
#     f_save = final_dir + '/sfigs/{}_{}.png'.format(sname, wstr)
#     plt.savefig(f_save, dpi=512)
#     plt.close()
#
# if use_mask:
#     lbl = 'rsq_crossterm_mask'
#     remove = pct_mask[0] * 100.0
# else:
#     lbl = 'rsq_crossterm'
#     remove = 0
# cross_term_correction_plot(varlist, rsqlist, rsdlist, pvar, lbl, remove)
#
#
# #%% SCATTER PLOT
# if use_mask:
#     lbl = 'ef_scatter_mask'
# else:
#     lbl = 'ef_scatter'
#
# for idx in [0, 1, 22, 3]:
#     plt.figure(figsize=(2.5, 2.5))
#     plt.subplots_adjust(top=0.96, right=0.96, bottom=0.2, left=0.22,
#                         hspace=0.5)
#
#     # set labels based on population index
#     pop = get_pop_list()[idx]
#     if 'mask' in lbl:
#         plab = '{} ({:.2f}% removed)'.format(pop, 100.0* pct_mask[idx])
#     else:
#         plab = pop
#
#     # select obs/pred
#     yi = nobs_list[idx]
#     fi = nprd_list[idx]
#     # get pred sorting index
#     fi_sidx = np.argsort(fi)
#     # sort data by pred
#     fi = fi[fi_sidx]
#     yi = yi[fi_sidx]
#     d_yf = yi-fi
#
#     # do linear regression on pred/obs
#     m, yint, rval, pval, stderr = linregress(fi, d_yf)
#     print pop, m
#
#     # do LOESS smoothing on pred/obs
#     wts = np.ones(shape=len(fi))
#     p_loess = predict_loess(fi, d_yf, wts, 0.1, fi)
#
#     # # calculate mean pred, residual in sliding windows
#     # # frac = 0.05
#     # # step = 0.5
#     # wsize = frac * (fi[-1] - fi[0])
#     # ssize = wsize * step
#     # mean_fi, mean_d_yf, num_sites = [], [], []
#     # for fi_i in np.arange(fi[0], fi[-1]+ssize-wsize, ssize):
#     #     fi_j = fi_i+wsize
#     #     istart = np.searchsorted(fi, fi_i)
#     #     iend = np.searchsorted(fi, fi_j)
#     #     if iend>istart:
#     #         mean_fi.append(fi[istart:iend].mean())
#     #         mean_d_yf.append(d_yf[istart:iend].mean())
#     #         num_sites.append(iend-istart)
#
#     # plt.subplot(5, 1, (2,5))
#     # plot scatter of pred vs. residuals
#     plt.plot(fi, d_yf,marker='o', ms=2.5,
#              markerfacecolor='None', markeredgecolor='k',
#              markeredgewidth=0.3, lw=0, alpha=0.75, label=plab)
#
#     # plot regression line
#     li = np.array([fi.min(), fi.max()])
#     plt.plot(li, li*m+yint, color='red', ls='--', lw=0.75,
#              label='linear regression')
#
#     # plot LOESS line
#     plt.plot(fi, p_loess, 'darkturquoise', lw=0.75, label='LOESS fit')
#
#     # # plot moving average line
#     # lbl += '_w{:.2f}_s{:.0f}'.format(wsize, 100*step)
#     # slbl = 'sliding window\nwin={:.2f}; slide={:.0f}%'.format(wsize, 100.0*step)
#     # plt.plot(mean_fi, mean_d_yf, 'darkturquoise', lw=0.75, marker='o',
#     #          label=slbl, ms=3)
#
#     plt.xlabel(r'$f_i$ ({} windows)'.format(wstr))
#     plt.ylabel(r'$y_i - f_i$ ({} windows)'.format(wstr))
#     plt.ylim(-1.4, 4)
#     plt.xlim(0.48, 2.15)
#     plt.legend(loc='upper left', frameon=1, framealpha=0.75,
#                facecolor='white', prop=dict(size=6))
#
#     # plt.subplot(5, 1, 1)
#     # plt.plot(mean_fi, num_sites, 'ko-', ms=2, label='sites/window',
#     #          lw=0.75)
#     # plt.xticks([])
#     # plt.xlim(0.48, 1.8)
#     # plt.legend(prop=dict(size=6))
#
#     f_save = final_dir + '/sfigs/{}_{}_{}.png'.format(pop, wstr, lbl)
#     plt.savefig(f_save, dpi=512)
#     plt.close()
#%% QUICK FUNCTION FOR LOADING SORT DATA
def get_sortpred_arrs(pop, nbins=100):
    # set path to results
    if pop == 'YRI':
        fldr = 'cadd94_gmask_mnb_378'
        fdir = final_dir + '/{}/'.format(fldr)
        sort_file = fdir + 'basic_sort_n{}.txt'.format(nbins)
    else:
        fldr = '{}_cadd94_gmask'.format(pop)
        fdir = final_dir + '/{}/'.format(fldr)
        sort_file = fdir + 'YRI_sorted.basic_sort_n{}.txt'.format(nbins)

    # load results data
    div, pi, pred = np.loadtxt(sort_file).T
    rst = cst_from_fldr(fldr)

    # normalize by pi0
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pi /= pi0
    pred /= pi0

    return pi, pred


#%% DIVERISTY DATA STRATIFIED BY YRI PREDICTIONS
def stratify_by_yri_map():
    # create new plot
    plt.figure(figsize=(3.25, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)
    # plot y=x line
    axmin, axmax = 0.55, 1.12
    plt.plot([axmin, axmax], [axmin, axmax], label=r'$y=x$',
             color='darkslategray', ls='--', alpha=0.65)
    # plot horizontal line at y=1
    plt.axhline(y=1, color='k', alpha=0.8, ls='-')
    plt.text(axmin + 0.01, 1.02, 'without linked selection', ha='left',
             va='center', fontsize=11)

    pops = 'YRI CEU GIH JPT MXL'.split()
    cols = 'darkorange steelblue darkturquoise purple fuchsia'.split()
    yripred = 0
    for i in xrange(5):
        pi, pred = get_sortpred_arrs(pops[i])
        if i == 0:
            yripred = pred
        # plot predicted vs. observed
        plt.plot(yripred, pi, marker='o', ms=5, markerfacecolor='None',
                 markeredgecolor=cols[i], markeredgewidth=0.9, lw=0,
                 alpha=0.75, label=pops[i])

    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
    ytick = np.arange(0.5, 1.2, 0.1)
    plt.yticks(ytick, x=0.02)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    xtick = np.arange(0.5, 1.2, 0.1)
    plt.xticks(xtick, y=0.02)
    plt.ylim(axmin, axmax)
    plt.xlim(axmin, axmax)
    plt.legend(loc='lower right', frameon=1, framealpha=0.75)

    f_save = final_dir + '/sfigs/compare.YRI.sorted.pops.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


stratify_by_yri_map()


#%%
def compared_stratified_b_bins():
    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)

    # YRI data
    yripi, yripred = get_sortpred_arrs('YRI')
    nyripi = yripi / yripi.mean()

    # remaining pop data to compare with YRI
    pops = 'CEU GIH JPT MXL'.split()
    cols = 'steelblue darkturquoise purple fuchsia'.split()
    for i in xrange(4):
        pi, pred = get_sortpred_arrs(pops[i])
        # plot predicted vs. observed
        lbl = '{} - YRI'.format(pops[i])
        diff = pi-yripi
        rat = pi/yripi
        npi = pi / np.mean(pi)
        ndiff = npi-nyripi
        nrat = npi/nyripi
        plt.plot(nrat, marker='o', ms=2.5, color=cols[i], lw=0.5,
                 # markeredgewidth=0.9,
                 alpha=0.75, label=lbl)

    plt.ylabel(r'$\mathrm{pop_{\pi/\pi_0} / YRI_{\pi/\pi_0}}$', labelpad=3)
    plt.xlabel(r'B bin', labelpad=3)
    plt.legend(loc='lower right', frameon=1, framealpha=0.75)

    f_save = final_dir + '/sfigs/pop.minus.YRI.across.Bbins.png'
    plt.savefig(f_save, dpi=512)
    plt.close()


compared_stratified_b_bins()


#%%


# yri = get_sortpred_arrs('YRI')
# ceu = get_sortpred_arrs('CEU')
# gih = get_sortpred_arrs('GIH')
# jpt = get_sortpred_arrs('JPT')
# mxl = get_sortpred_arrs('MXL')
#
# all_arr = [yri, ceu, gih, jpt, mxl]


#%%
def compare_population_diversity_across_b_bins(nbins=100, normed=False,
                                               use_bvals=False, span=0.1):
    # create new plot
    plt.figure(figsize=(6.5, 3.25))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.15, bottom=0.12)
    ax1 = plt.subplot(111)
    format_panels(ax1)
    plt.title('LOESS span = {:.0f}%'.format(100*span), y=0.75)

    # YRI data
    yripi, yripred = get_sortpred_arrs('YRI', nbins)
    # get normalized prediction and pi values
    nyripi = yripi / yripi.mean()
    # nyripred = yripred / yripred.mean()

    # remaining pop data to compare with YRI
    pops = 'CEU GIH JPT MXL'.split()
    cols = 'steelblue darkturquoise purple fuchsia'.split()
    loess_points = []
    for i in xrange(4):
        pi, pred = get_sortpred_arrs(pops[i], nbins)
        # lbl = '{} (LOESS span = 10%)'.format(pops[i])
        # diff = pi-yripi
        # rat = pi/yripi
        # normalize pi by its mean
        npi = pi / np.mean(pi)
        # take absolute difference between YRI and other population across bins
        ndiff = npi-nyripi
        # scale the difference by normalized B in each bin if normed flag True
        if normed:
            ndiff /= yripred
        # get LOESS fit for the points
        wts = np.ones(shape=len(ndiff))
        ndiff_loess = predict_loess(yripred, ndiff, wts, span, yripred)
        loess_points.append(ndiff_loess)
        if use_bvals:
            plt.plot(yripred, ndiff, marker='o', ms=1.5, color=cols[i], lw=0,
                     alpha=0.5)
        else:
            plt.plot(ndiff, marker='o', ms=1.5, color=cols[i], lw=0, alpha=0.5)

    # plot lines on top of all points
    for i in xrange(4):
        ndiff_loess = loess_points[i]
        if use_bvals:
            plt.plot(yripred, ndiff_loess, color=cols[i], label=pops[i],
                     lw=1.5)
        else:
            plt.plot(ndiff_loess, color=cols[i], label=pops[i], lw=1.5)

    ylbl = r'$\mathrm{pop_{\pi/\pi_0} - YRI_{\pi/\pi_0}}$'
    if normed:
        ylbl += ' (normalized by B)'
    plt.ylabel(ylbl, labelpad=3)
    xlbl = 'B value' if use_bvals else 'B bin'
    plt.xlabel(xlbl, labelpad=3)
    plt.legend(loc='lower right', frameon=1, framealpha=0.75,
               facecolor='white')
    plt.ylim(-0.125, 0.125)
    # if normed:
    #     plt.ylim(-0.125, 0.125)
    # else:
    #     plt.ylim(-0.4, 0.3)
    pltlbl = 'compare_population_diversity_across_b_bins'
    pltlbl += '_span_{:.0f}pct'.format(100*span)
    if normed:
        pltlbl += '_normalized'
    if use_bvals:
        pltlbl += '_use_bvals'
    f_save = final_dir + '/sfigs/{}.png'.format(pltlbl)

    plt.savefig(f_save, dpi=512)
    plt.close()


# compare_population_diversity_across_b_bins(1000, span=0.1)
# compare_population_diversity_across_b_bins(1000, normed=True)
# compare_population_diversity_across_b_bins(1000, use_bvals=True)
# compare_population_diversity_across_b_bins(1000, normed=True, use_bvals=True)
compare_population_diversity_across_b_bins(1000, span=0.2)


#%%
