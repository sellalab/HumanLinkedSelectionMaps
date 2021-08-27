__author__ = 'davidmurphy'



import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from classes.runstruct import ChromStruct, root_dir
from data_processing.data_tools import randint_unique
from data_processing.functions import rsquared_function
from figures.common_functions import format_panels, get_pop_array, \
    get_pop_array_2


final_dir = root_dir + '/result/final_files'


#%% OVERFITTING RELATIVE DIFFERENCES PLOT
def overfitting_plot(an1, an2):
    f_1 = final_dir + '/{an}/rsq.log'.format(an=an1)
    f_jk1 = final_dir + '/{an}_jackknife_results/{an}.jkrsq.log'.format(an=an1)
    f_2 = final_dir + '/{an}/rsq.log'.format(an=an2)
    f_jk2 = final_dir + '/{an}_jackknife_results/{an}.jkrsq.log'.format(an=an2)

    groups = [[f_1, f_jk1], [f_2, f_jk2]]
    lbls = ['phastCons', 'CADD']
    ttls = ['phastCons 6%', 'CADD 6%']
    cols = ['darkorange', 'dodgerblue']
    lets = ['a', 'b']
    lpos = [0.105, 0.612]
    plt.figure(figsize=(5.4, 2.16))
    plt.subplots_adjust(bottom=0.17, top=0.9, right=0.995, left=0.1,
                        wspace=0.3)
    for i in xrange(2):
        ax1 = plt.subplot(1, 2, i+1)
        format_panels(ax1)
        plt.title(ttls[i], loc='center')
        f, fjk = groups[i]
        w, r = np.loadtxt(f)[:16].T
        jw, jr = np.loadtxt(fjk)[:16].T
        # print zip(r, jr)
        d_rel = 100.0 * (r-jr) / r
        # print d_rel
        plt.plot(np.log10(w), d_rel, marker='o', lw=0, color=cols[i],
                 label=lbls[i])
        plt.xlabel('window size (log-scale)', labelpad=2)
        xtck = [4, 4.5, 5, 5.5, 6]
        xstr = [r'$10^{%.1f}$' % x if not x % 1 else '' for x in xtck]
        plt.xticks(xtck, xstr, y=0.04)
        # plt.ylim(0.05, 0.65)
        if i == 0:
            plt.ylabel(r'relative difference (%)', labelpad=2)
        plt.yticks(x=0.02)

        plt.text(lpos[i], 0.83, lets[i], transform=plt.gcf().transFigure,
                 fontweight='bold')

        # plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], x=0.02)
        # plt.legend(loc='center left', borderpad=-0.1, handletextpad=-0.1,
        #            borderaxespad=-0.1)

    fsave = final_dir + '/sfigs/overfit.relative.diff.png'
    plt.savefig(fsave, dpi=512)
    plt.close()


overfitting_plot('fish_cons94_new', 'cadd94_gmask_v1.6_without_bstat')


#%% JACKKNIFE PLOT
def jackknife_variance(jk_estimates):
    n = len(jk_estimates)
    jk_mean = np.average(jk_estimates, axis=0)
    jk_var = np.sum(np.power(jk_estimates-jk_mean, 2), axis=0) * (n-1) / n

    return jk_var


def jackknife_stats():
    """get parameters with CI from jackknife results"""
    # labels and colors for two main results
    colors = ['darkorange', 'dodgerblue']
    labels = ['phastCons', 'CADD']
    u0 = 1.4

    # get jackknife results from each main run
    results = []
    for anno in ['fish_cons94_new', 'cadd94_gmask_v1.6_without_bstat']:
        save_dir = final_dir + '/{an}_jackknife_results/'.format(an=anno)
        # create list of filenames for saving or loading presaved data
        ftokens = ['pmf.npy', 'udl.npy', 'pi0.npy', 'clh.npy']
        # load data arrays
        dtfiles = [save_dir + ft for ft in ftokens]
        # pmf, udl, pi0, clh = [np.load(f) for f in dtfiles]
        results.append([np.load(f) for f in dtfiles])

    # create the plot
    plt.figure(figsize=(6.5, 2.16))
    # plt.suptitle('jackknife 2Mb windows')
    plt.subplots_adjust(left=0.07, top=0.97, wspace=1.2, right=0.995,
                        bottom=0.17)

    # plot DFE
    ax1 = plt.subplot(1, 5, (1,3))
    format_panels(ax1)
    xax = np.arange(6)
    width = 0.4
    offset = -0.2
    for (i, res) in enumerate(results):
        pmf = res[0] / u0
        # average the jackknife PMFs
        pmf_avg = np.average(pmf, axis=0)
        # get the jackknife variance for error bars
        pmf_erange = np.sqrt(jackknife_variance(pmf)) / 2.0
        plt.bar(xax+offset, pmf_avg, width=width, yerr=pmf_erange,
                ecolor='k', color=colors[i], label=labels[i])
        offset += width
    xtck = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0]
    plt.xticks(xax, [r'$10^{%.1f}$' % x for x in xtck], y=0.03)
    plt.xlabel('deleterious fitness effect', labelpad=2)
    plt.yticks(x=0.02)
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=2)
    plt.legend(loc='upper left', frameon=1, framealpha=0.75, facecolor='white',
               handlelength=1, labelspacing=0.25)
    plt.ylim(0, 0.45)
    # make small ticks between selection coefs
    for xx in xax[:-1]:
        plt.axvline(xx+0.5, 0, 0.2, lw=1, ls='--', color='k')

    # plot udel
    ax2 = plt.subplot(154)
    format_panels(ax2)
    for (i, res) in enumerate(results):
        udl = res[1] / u0
        udl_erange = np.sqrt(jackknife_variance(udl)) / 2.0
        plt.bar(i, udl.mean(), width=1, yerr=udl_erange, ecolor='k',
                color=colors[i])
    plt.ylabel(r'mutation rate (units of $u_0$)', labelpad=3)
    plt.xticks([])

    # plot pi0
    ax3 = plt.subplot(155)
    format_panels(ax3)
    for (i, res) in enumerate(results):
        pi0 = res[2]
        pi0_erange = np.sqrt(jackknife_variance(pi0)) / 2.0
        plt.bar(i, pi0.mean(), width=1, yerr=pi0_erange, ecolor='k',
                color=colors[i])
    plt.ylabel(r'diversity reduction ($\bar{\pi}/\pi_0)$', labelpad=1)
    plt.xticks([])
    # plt.ylim(pi0.min()-2*pi0_erange, pi0.max()+2*pi0_erange)

    # save high res figure
    fsave = final_dir + '/sfigs/jackknife.params.png'
    plt.savefig(fsave, dpi=512)
    plt.close()


# anno = 'fish_cons94_gmask'
jackknife_stats()


#%% FIXED RESAMPLING STRATEGY (VIA GUY AMSTER) FOR NULL DISTRIBUTION
def create_null_dist(win, n_resamples):
    """
    create a null distribution of the difference between R^2 values
    by randomly constructing maps from a shuffle of map A and map B
    """
    # get both sets of results
    obs1, prd1 = get_pop_array(win, 'YRI', 'fish_cons94_new')
    obs2, prd2 = get_pop_array(win, 'YRI', 'cadd94_gmask_v1.6_without_bstat')

    # set size of random integer space to sample data indices from
    imax = len(obs1)

    # PATCH: randomly trim windows to match sizes
    if len(obs1) > len(obs2):
        imax = len(obs2)
        size_diff = len(obs1) - len(obs2)
        msk = np.ones(shape=len(obs1), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs1))
        msk[rnd_rm] = False
        rsqbefore = rsquared_function(obs1, prd1)
        obs1, prd1 = obs1[msk], prd1[msk]
        rsqafter  = rsquared_function(obs1, prd1)
        msg = 'size_diff={} rsq_before={} rsq_after={}'
        print msg.format(size_diff, rsqbefore, rsqafter)

    if len(obs2) > len(obs1):
        size_diff = len(obs2) - len(obs1)
        msk = np.ones(shape=len(obs2), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs2))
        msk[rnd_rm] = False
        rsqbefore = rsquared_function(obs2, prd2)
        obs2, prd2 = obs2[msk], prd2[msk]
        rsqafter  = rsquared_function(obs2, prd2)
        msg = 'size_diff={} rsq_before={} rsq_after={}'
        print msg.format(size_diff, rsqbefore, rsqafter)

    assert len(obs1) == len(obs2)

    # get the true R^2 for each run
    true_1 = rsquared_function(obs1, prd1)
    true_2 = rsquared_function(obs2, prd2)

    # lists for R^2 results from bootstraps
    rsq_1, rsq_2 = [], []

    # group smaller windows into approximate Mb windows
    if win != 1e6:
        # determine the number of groups to collect data into
        n_groups = int(1e6/win)
        # determine the length of windows to include in groups
        n_window = int(imax/n_groups)
        # put grouped windows into a new array
        group_obs1 = np.zeros(shape=(n_window, n_groups))
        group_obs2 = np.zeros(shape=(n_window, n_groups))
        group_prd1 = np.zeros(shape=(n_window, n_groups))
        group_prd2 = np.zeros(shape=(n_window, n_groups))
        for i in xrange(n_window):
            # j = i+1
            istart = i*n_groups
            iend = istart + n_groups
            group_obs1[i] = obs1[istart:iend]
            group_obs2[i] = obs2[istart:iend]
            group_prd1[i] = prd1[istart:iend]
            group_prd2[i] = prd2[istart:iend]

        for _ in xrange(n_resamples):
            # create random T/F array to select from each map
            msk = np.random.choice([True, False], size=n_window, p=[0.5, 0.5])

            # for obs and pred A, take T for map1 and F for map2
            i1, i2 = np.where(msk)[0], np.where(~msk)[0]
            o1 = np.concatenate((np.concatenate(group_obs1[i1]),
                                 np.concatenate(group_obs2[i2])))
            p1 = np.concatenate((np.concatenate(group_prd1[i1]),
                                 np.concatenate(group_prd2[i2])))

            # idx = np.random.randint(0, n_window, size=n_window)
            # o1 = np.concatenate(group_obs1[idx])
            # p1 = np.concatenate(group_prd1[idx])
            r1 = rsquared_function(o1, p1)
            rsq_1.append(r1)

            # for obs and pred B, take the opposite set of windows
            o2 = np.concatenate((np.concatenate(group_obs1[i2]),
                                 np.concatenate(group_obs2[i1])))
            p2 = np.concatenate((np.concatenate(group_prd1[i2]),
                                 np.concatenate(group_prd2[i1])))
            r2 = rsquared_function(o2, p2)
            rsq_2.append(r2)
    else:
        # create N randomly shuffles of maps A and B
        for _ in xrange(n_resamples):
            # create random T/F array to select from each map
            msk = np.random.choice([True, False], size=imax, p=[0.5, 0.5])

            # for obs and pred A, take the T for obs1 and the F for obs2
            obs_a = np.concatenate((obs1[msk], obs2[~msk]))
            prd_a = np.concatenate((prd1[msk], prd2[~msk]))

            # for obs and pred B, take the opposite set of windows
            obs_b = np.concatenate((obs1[~msk], obs2[msk]))
            prd_b = np.concatenate((prd1[~msk], prd2[msk]))

            # calculate R^2 for map A and map B and take the difference
            rsq_a = rsquared_function(obs_a, prd_a)
            rsq_1.append(rsq_a)
            rsq_b = rsquared_function(obs_b, prd_b)
            rsq_2.append(rsq_b)

    return true_1, true_2, np.array(rsq_1), np.array(rsq_2)


# null1mb = create_null_dist(1e6, 10000)
# tfish, tcad, rfish, rcad = resample_maps(1e6, 10000)
null_dists = [create_null_dist(w, 10000) for w in [1.1e4, 1.25e5, 1e6]]



#%% NEW VERSION USING GUY A'S METHOD
def plot_rsq_bootstraps_3(null_dists, nbins=100):
    """plot the distribution of boostrap R^2 values and their differences"""
    stitles = ['10 kb windows', '100 kb windows', '1 Mb windows']
    wins = ['10kb', '100kb', '1Mb']
    colors = ['darkorange', 'dodgerblue']
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(left=0.085, top=0.9, right=0.99, wspace=0.3,
                        bottom=0.19, hspace=0.4)

    # plot each set of null distributions with the real difference indicated
    for i in xrange(3):
        ax = plt.subplot(1, 3, i+1)
        plt.title(stitles[i], loc='center', y=0.975)
        format_panels(ax)

        # get data from null dist list
        tfish, tcad, rfish, rcad = null_dists[i]

        # get the histogram of R^2 differences
        true_diff = tcad-tfish
        rdiff = rcad-rfish
        cnts, bins = np.histogram(rdiff, bins=nbins)
        cnts = cnts.astype(float) / 1e4

        # set params for bar plot of differences
        widths = bins[1:] - bins[:-1]
        xi = bins[:-1] + widths/2.0
        mneg = (xi<true_diff)
        mpos = (xi>=true_diff)
        plt.bar(xi[mneg], cnts[mneg], widths[mneg], color='gray', alpha=0.75)
        plt.bar(xi[mpos], cnts[mpos], widths[mpos], color='red', alpha=0.75)
        # plt.axvline(0, color='k', ls=':', lw=0.5)

        # calculate the p-value based on the values in the distribution
        r_plus1 = 1.0 * np.sum(rdiff >= true_diff) + 1
        n_plus1 = len(rcad) + 1
        pval = r_plus1 / n_plus1

        pstr = r'$p={:.3f}$'.format(pval)
        plt.text(0, 0.035, pstr, fontsize=10, ha='center', va='bottom')

        ymax = 0.04
        ytcks = np.array([100, 200, 300, 400]) / 1e4
        plt.yticks(ytcks, x=0.04)
        plt.ylim(0, ymax)
        if i == 0:
            plt.ylabel(r'relative frequency', labelpad=2)
        # plt.ylim(0, ymax)
        plt.xticks(y=0.04)
        plt.xlabel(r'$\mathrm{\Delta R^2}$', labelpad=2)

    # plt.legend()
    fsave = final_dir + '/sfigs/CADD-phastCons-permutation.rsq.png'
    plt.savefig(fsave, dpi=512)
    plt.close()


plot_rsq_bootstraps_3(null_dists)


#%% RESAMPLE, COMPARE R^2 AT 1MB FOR DIFFERENT DEPTHS/PERCENTS (GUY A. METHOD)
def resample_maps_3(fldr1, fldr2, win, n_resamples):
    """resample from main best models and calculate R^2"""
    # get both sets of results
    obs1, prd1 = get_pop_array_2(win, fldr1, 'self')
    obs2, prd2 = get_pop_array_2(win, fldr2, 'self')

    # set size of random integer space to sample data indices from
    imax = len(obs1)

    # PATCH: randomly trim windows to match sizes
    if len(obs1) > len(obs2):
        imax = len(obs2)
        size_diff = len(obs1) - len(obs2)
        msk = np.ones(shape=len(obs1), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs1))
        msk[rnd_rm] = False
        rsqbefore = rsquared_function(obs1, prd1)
        obs1, prd1 = obs1[msk], prd1[msk]
        rsqafter = rsquared_function(obs1, prd1)
        msg = 'size_diff={} rsq_before={} rsq_after={}'
        print msg.format(size_diff, rsqbefore, rsqafter)

    if len(obs2) > len(obs1):
        size_diff = len(obs2) - len(obs1)
        msk = np.ones(shape=len(obs2), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs2))
        msk[rnd_rm] = False
        rsqbefore = rsquared_function(obs2, prd2)
        obs2, prd2 = obs2[msk], prd2[msk]
        rsqafter = rsquared_function(obs2, prd2)
        msg = 'size_diff={} rsq_before={} rsq_after={}'
        print msg.format(size_diff, rsqbefore, rsqafter)

    assert len(obs1) == len(obs2)

    # get the true R^2 for each run
    true_1 = rsquared_function(obs1, prd1)
    true_2 = rsquared_function(obs2, prd2)

    # lists for R^2 results from bootstraps
    rsq_1, rsq_2 = [], []

    # group smaller windows into approximate Mb windows
    if win != 1e6:
        # determine the number of groups to collect data into
        n_groups = int(1e6 / win)
        # determine the length of windows to include in groups
        n_window = int(imax / n_groups)
        # put grouped windows into a new array
        group_obs1 = np.zeros(shape=(n_window, n_groups))
        group_obs2 = np.zeros(shape=(n_window, n_groups))
        group_prd1 = np.zeros(shape=(n_window, n_groups))
        group_prd2 = np.zeros(shape=(n_window, n_groups))
        for i in xrange(n_window):
            # j = i+1
            istart = i * n_groups
            iend = istart + n_groups
            group_obs1[i] = obs1[istart:iend]
            group_obs2[i] = obs2[istart:iend]
            group_prd1[i] = prd1[istart:iend]
            group_prd2[i] = prd2[istart:iend]

        for _ in xrange(n_resamples):
            # create random T/F array to select from each map
            msk = np.random.choice([True, False], size=n_window, p=[0.5, 0.5])

            # for obs and pred A, take T for map1 and F for map2
            i1, i2 = np.where(msk)[0], np.where(~msk)[0]
            o1 = np.concatenate((np.concatenate(group_obs1[i1]),
                                 np.concatenate(group_obs2[i2])))
            p1 = np.concatenate((np.concatenate(group_prd1[i1]),
                                 np.concatenate(group_prd2[i2])))

            # idx = np.random.randint(0, n_window, size=n_window)
            # o1 = np.concatenate(group_obs1[idx])
            # p1 = np.concatenate(group_prd1[idx])
            r1 = rsquared_function(o1, p1)
            rsq_1.append(r1)

            # for obs and pred B, take the opposite set of windows
            o2 = np.concatenate((np.concatenate(group_obs1[i2]),
                                 np.concatenate(group_obs2[i1])))
            p2 = np.concatenate((np.concatenate(group_prd1[i2]),
                                 np.concatenate(group_prd2[i1])))
            r2 = rsquared_function(o2, p2)
            rsq_2.append(r2)
    else:
        # create N randomly shuffles of maps A and B
        for _ in xrange(n_resamples):
            # create random T/F array to select from each map
            msk = np.random.choice([True, False], size=imax, p=[0.5, 0.5])

            # for obs and pred A, take the T for obs1 and the F for obs2
            obs_a = np.concatenate((obs1[msk], obs2[~msk]))
            prd_a = np.concatenate((prd1[msk], prd2[~msk]))

            # for obs and pred B, take the opposite set of windows
            obs_b = np.concatenate((obs1[~msk], obs2[msk]))
            prd_b = np.concatenate((prd1[~msk], prd2[msk]))

            # calculate R^2 for map A and map B and take the difference
            rsq_a = rsquared_function(obs_a, prd_a)
            rsq_1.append(rsq_a)
            rsq_b = rsquared_function(obs_b, prd_b)
            rsq_2.append(rsq_b)

    return true_1, true_2, np.array(rsq_1), np.array(rsq_2)

    # # set size of random integer space to sample data indices from
    # imax = len(obs1)
    #
    # # PATCH
    # if len(obs1) > len(obs2):
    #     imax = len(obs2)
    #     size_diff = len(obs1) - len(obs2)
    #     msk = np.ones(shape=len(obs1), dtype=bool)
    #     rnd_rm = randint_unique(size_diff, len(obs1))
    #     msk[rnd_rm] = False
    #     print size_diff, rsquared_function(obs1, prd1)
    #     obs1, prd1 = obs1[msk], prd1[msk]
    #     print rsquared_function(obs1, prd1)
    # if len(obs2) > len(obs1):
    #     size_diff = len(obs2) - len(obs1)
    #     msk = np.ones(shape=len(obs2), dtype=bool)
    #     rnd_rm = randint_unique(size_diff, len(obs2))
    #     msk[rnd_rm] = False
    #     print size_diff, rsquared_function(obs2, prd2)
    #     obs2, prd2 = obs2[msk], prd2[msk]
    #     print rsquared_function(obs2, prd2)
    #
    # assert len(obs1) == len(obs2)
    #
    # # get the true R^2 for each run
    # true_1 = rsquared_function(obs1, prd1)
    # true_2 = rsquared_function(obs2, prd2)
    #
    # # lists for R^2 results from bootstraps
    # rsq_1, rsq_2 = [], []
    #
    # # group smaller windows into approximate Mb windows
    # if win != 1e6:
    #     # determine the number of groups to collect data into
    #     n_groups = int(1e6/win)
    #     # determine the length of windows to include in groups
    #     n_window = int(imax/n_groups)
    #     # put grouped windows into a new array
    #     group_obs1 = np.zeros(shape=(n_window, n_groups))
    #     group_obs2 = np.zeros(shape=(n_window, n_groups))
    #     group_prd1 = np.zeros(shape=(n_window, n_groups))
    #     group_prd2 = np.zeros(shape=(n_window, n_groups))
    #     for i in xrange(n_window):
    #         # j = i+1
    #         istart = i*n_groups
    #         iend = istart + n_groups
    #         group_obs1[i] = obs1[istart:iend]
    #         group_obs2[i] = obs2[istart:iend]
    #         group_prd1[i] = prd1[istart:iend]
    #         group_prd2[i] = prd2[istart:iend]
    #
    #     for _ in xrange(n_resamples):
    #         idx = np.random.randint(0, n_window, size=n_window)
    #         o1 = np.concatenate(group_obs1[idx])
    #         p1 = np.concatenate(group_prd1[idx])
    #         r1 = rsquared_function(o1, p1)
    #         rsq_1.append(r1)
    #
    #         o2 = np.concatenate(group_obs2[idx])
    #         p2 = np.concatenate(group_prd2[idx])
    #         r2 = rsquared_function(o2, p2)
    #         rsq_2.append(r2)
    # else:
    #     # do N resamples of each data set anc calculate the R^2 values
    #     for _ in xrange(n_resamples):
    #         # sample indices resample with replacement
    #         idx = np.random.randint(0, imax, size=imax)
    #         # calculate R^2 from each set of data using the same idx
    #         r1 = rsquared_function(obs1[idx], prd1[idx])
    #         rsq_1.append(r1)
    #         r2 = rsquared_function(obs2[idx], prd2[idx])
    #         rsq_2.append(r2)
    #
    # return true_1, true_2, np.array(rsq_1), np.array(rsq_2)


#%% GET DATA FOR R^2 COMPARISONS (GUY A. METHOD)
fl1 = 'fish_cons94_new'
spec = 'ape primate prosimian euarchontoglires laurasiatheria mammal'.split()
n_boots = 10000
spec_pvals = []
print 'species'
for sp in spec:
    if sp == 'ape':
        fl2 = '{}_cons94_new'.format(sp)
    else:
        fl2 = '{}_cons94_new'.format(sp)
    t1, t2, bs1, bs2 = resample_maps_3(fl1, fl2, 1e6, n_boots)
    pval = 1.0 * (np.sum(bs1-bs2>=t1-t2)+1.0) / (n_boots+1.0)
    print sp, t1, t2, pval
    spec_pvals.append(pval)
print '---'
#%%
fl1 = 'cadd94_gmask_v1.6_without_bstat'
pct = [91, 92, 93, 95, 96, 97, 98]
n_boots = 10000
cadd_pvals = []
print 'CADD'
for p in pct:
    fl2 = 'cadd{}_gmask_v1.6_without_bstat'.format(p)
    t1, t2, bs1, bs2 = resample_maps_3(fl1, fl2, 1e6, n_boots)
    pval = 1.0 * (np.sum(bs1-bs2>=t1-t2)+1.0) / (n_boots+1.0)
    # pval = 1.0 * np.sum(bs1<=bs2) / n_boots
    print p, t1, t2, pval
    cadd_pvals.append(pval)
print '---'
#%%
fl1 = 'fish_cons94_new'
pct = [91, 92, 93, 95, 96, 97, 98]
n_boots = 10000
fish_pvals = []
print 'fish percent'
for p in pct:
    fl2 = 'fish_cons{}_new'.format(p)
    t1, t2, bs1, bs2 = resample_maps_3(fl1, fl2, 1e6, n_boots)
    pval = 1.0 * (np.sum(bs1-bs2>=t1-t2)+1.0) / (n_boots+1.0)
    # pval = 1.0 * np.sum(bs1<=bs2) / n_boots
    print p, t1, t2, pval
    fish_pvals.append(pval)


#%% PLOT PVAL FOR R^2 COMPARISONS (GUY A. METHOD)
def plot_pvals(pvals, labels, xlabs):
    plt.figure(figsize=(6.5, 2.16))
    plt.subplots_adjust(top=0.995, left=0.08, right=0.995, bottom=0.3,
                        wspace=0.15)
    txt_x = [0.083, 0.404, 0.72]
    letters = ['a', 'b', 'c']
    for i in xrange(len(pvals)):
        ax = plt.subplot(1, 3, i+1)
        n_stars = 0
        format_panels(ax)
        pvs, lbs = pvals[i], labels[i]
        xi = np.arange(len(pvs))
        plt.axhline(-np.log10(0.05), ls='--', color='k')
        for j in xrange(len(xi)):
            if pvs[j] > 0:
                plt.bar(xi[j], -np.log10(pvs[j]), color='gray', alpha=0.75)
            else:
                plt.bar(xi[j], 5, color='gray', alpha=0.75)
                #
                # if n_stars:
                #     plt.plot(xi[j], 3, '*', color='gray')
                # else:
                #     plt.plot(xi[j], 3, '*', label='*p = 0', color='gray')
                # n_stars += 1
        plt.xticks(xi, lbs, rotation=45, y=0.04)
        plt.xlabel(xlabs[i], labelpad=2)
        if i == 0:
            plt.ylabel('significance level\n' +
                       r'$\mathrm{(-log_{10}}$(p-value))', labelpad=2)
        y_tck = [0, 1, 2, 3, 4]
        plt.yticks(y_tck, x=0.04)
        plt.ylim(0, 4.2)
        plt.text(txt_x[i], 0.92, letters[i], transform=plt.gcf().transFigure,
                 fontweight='bold')
        if i == 0:
            plt.text(-0.4, 1.5, 'p=0.05', va='bottom', ha='left')
            # bbox = dict(facecolor='white', alpha=0.5, lw=0))
            # plt.legend(loc='center left', frameon=1, framealpha=0.5,
        #            facecolor='white')

    f_save = final_dir + '/sfigs/updateMarch2021-permutation-pvalue_distributions.png'
    plt.savefig(f_save, dpi=512)


pct_labls = [2, 3, 4, 5, 7, 8, 9]
pvals = [spec_pvals, fish_pvals[::-1], cadd_pvals[::-1]]
labels = [[s[:4] for s in spec], pct_labls, pct_labls]
xlabs = ['phylogenetic depth', '% conserved', '% CADD']
plot_pvals(pvals, labels, xlabs)


#%% RESAMPLE FOR R^2 BOOTSTRAP (OLD - NOT USED)
def resample_maps(win, n_resamples):
    """resample from main best models and calculate R^2"""
    # get both sets of results
    obs1, prd1 = get_pop_array(win, 'YRI', 'fish_cons94_gmask')
    obs2, prd2 = get_pop_array(win, 'YRI', 'cadd94_gmask')

    # set size of random integer space to sample data indices from
    imax = len(obs1)

    # PATCH
    if len(obs1) > len(obs2):
        imax = len(obs2)
        size_diff = len(obs1) - len(obs2)
        msk = np.ones(shape=len(obs1), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs1))
        msk[rnd_rm] = False
        print size_diff, rsquared_function(obs1, prd1)
        obs1, prd1 = obs1[msk], prd1[msk]
        print rsquared_function(obs1, prd1)
    if len(obs2) > len(obs1):
        size_diff = len(obs2) - len(obs1)
        msk = np.ones(shape=len(obs2), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs2))
        msk[rnd_rm] = False
        print size_diff, rsquared_function(obs2, prd2)
        obs2, prd2 = obs2[msk], prd2[msk]
        print rsquared_function(obs2, prd2)

    assert len(obs1) == len(obs2)

    # get the true R^2 for each run
    true_1 = rsquared_function(obs1, prd1)
    true_2 = rsquared_function(obs2, prd2)

    # lists for R^2 results from bootstraps
    rsq_1, rsq_2 = [], []

    # group smaller windows into approximate Mb windows
    if win != 1e6:
        # determine the number of groups to collect data into
        n_groups = int(1e6/win)
        # determine the length of windows to include in groups
        n_window = int(imax/n_groups)
        # put grouped windows into a new array
        group_obs1 = np.zeros(shape=(n_window, n_groups))
        group_obs2 = np.zeros(shape=(n_window, n_groups))
        group_prd1 = np.zeros(shape=(n_window, n_groups))
        group_prd2 = np.zeros(shape=(n_window, n_groups))
        for i in xrange(n_window):
            # j = i+1
            istart = i*n_groups
            iend = istart + n_groups
            group_obs1[i] = obs1[istart:iend]
            group_obs2[i] = obs2[istart:iend]
            group_prd1[i] = prd1[istart:iend]
            group_prd2[i] = prd2[istart:iend]

        for _ in xrange(n_resamples):
            idx = np.random.randint(0, n_window, size=n_window)
            o1 = np.concatenate(group_obs1[idx])
            p1 = np.concatenate(group_prd1[idx])
            r1 = rsquared_function(o1, p1)
            rsq_1.append(r1)

            o2 = np.concatenate(group_obs2[idx])
            p2 = np.concatenate(group_prd2[idx])
            r2 = rsquared_function(o2, p2)
            rsq_2.append(r2)
    else:
        # do N resamples of each data set anc calculate the R^2 values
        for _ in xrange(n_resamples):
            # sample indices resample with replacement
            idx = np.random.randint(0, imax, size=imax)
            # calculate R^2 from each set of data using the same idx
            r1 = rsquared_function(obs1[idx], prd1[idx])
            rsq_1.append(r1)
            r2 = rsquared_function(obs2[idx], prd2[idx])
            rsq_2.append(r2)

    return true_1, true_2, np.array(rsq_1), np.array(rsq_2)


# tfish, tcad, rfish, rcad = resample_maps(1e6, 10000)
bootstraps = [resample_maps(w, 10000) for w in [1.1e4, 1.25e5, 1e6]]


#%% BOOTSTRAP PLOTS COMPARING R^2 (FIRST VERSION -- NOT USED IN S7)
def plot_rsq_bootstraps(rsq_list, nbins=100):
    """plot the distribution of boostrap R^2 values and their differences"""
    labels = ['99-vert', 'CADD']
    colors = ['darkorange', 'dodgerblue']

    plt.figure(figsize=(6.5, 4.32))
    plt.subplots_adjust(left=0.1, top=0.995, right=0.995, wspace=0.4,
                        bottom=0.1, hspace=0.4)
    for (i, rsq_pair) in enumerate(rsq_list):
        # plot histogram of R^2 values
        ax = plt.subplot(2, 4, 2*i+1)
        format_panels(ax)
        tfish, tcad, rfish, rcad = rsq_pair
        # plot distribution of fish and cadd results
        # cnts, bins = np.histogram(rfish, bins=nbins)
        # cnts = 1.0 * cnts / np.sum(cnts)
        # widths = bins[1:] - bins[:-1]
        # xi = bins[:-1] + widths / 2.0
        plt.hist(rfish, bins=nbins, histtype='step', color=colors[0],
                 label=labels[0], alpha=0.75, lw=0.75)
        # plt.step(xi, cnts, color=colors[0], label=labels[0], alpha=0.75,
        #          lw=0.75)
        # cnts, bins = np.histogram(rcad, bins=nbins)
        # cnts = 1.0 * cnts / np.sum(cnts)
        # widths = bins[1:] - bins[:-1]
        # xi = bins[:-1] + widths / 2.0
        # plt.step(xi, cnts, color=colors[1], label=labels[1], alpha=0.75,
        #          lw=0.75)
        plt.hist(rcad, bins=nbins, histtype='step', color=colors[1],
                 label=labels[1], alpha=0.75, lw=0.75)
        # plot true values from each
        plt.axvline(x=tfish, color=colors[0], ls='--', lw=1)
        plt.axvline(x=tcad, color=colors[1], ls='--', lw=1)
        # plt.xticks(rotation=90)
        # plot histogram of differences in R^2
        ax = plt.subplot(2, 4, 2*i+2)
        format_panels(ax)
        cnts, bins = np.histogram(rfish-rcad, bins=nbins)
        widths = bins[1:] - bins[:-1]
        xi = bins[:-1] + widths/2.0
        mneg = (xi<0)
        mpos = (xi>=0)
        plt.bar(xi[mneg], cnts[mneg], widths[mneg], color=colors[1])
        plt.bar(xi[mpos], cnts[mpos], widths[mpos], color=colors[0])
        plt.axvline(0, color='k', ls=':', lw=0.5)
        # plt.hist(rfish-rcad, bins=100, color='gray', label='99-vert - CADD')
        # plt.xticks(rotation=90)
        # plt.yticks(color='none')
        print 1.0*np.sum(rfish-rcad>0)/len(rcad)

    fsave = final_dir + '/sfigs/bootstrap.rsq.png'
    plt.savefig(fsave, dpi=512)
    plt.close()


#%% BOOTSTRAP PLOTS COMPARING R^2 (OLD - NOT USED)
def plot_rsq_bootstraps_2(rsq_pair, idx, nbins=100):
    """plot the distribution of boostrap R^2 values and their differences"""
    stitles = ['A) 10kb windows', 'B) 100kb windows', 'C) 1Mb windows']
    wins = ['10kb', '100kb', '1Mb']
    colors = ['darkorange', 'dodgerblue']
    xt1 = [[0.145, 0.15, 0.155, 0.16], [0.38, 0.4, 0.42], [0.55, 0.6, 0.65]]
    xl1 = [(0.1455, 0.1595), (0.365, 0.425), (0.54, 0.66)]
    xt2 = [[-0.008, 0, 0.008], [-0.02, 0, 0.02], [-0.01, 0, 0.01]]
    xl2 = [(-0.011, 0.011), (-0.03, 0.03), (-0.016, 0.016)]
    tpos = [5e-5, 5e-5, 5e-5]
    tplus = [0.00035, 0.0015, 0.003]
    ymax = 375
    plt.figure(figsize=(5.4, 2.16))
    plt.subplots_adjust(left=0.085, top=0.9, right=0.99, wspace=0.3,
                        bottom=0.19, hspace=0.4)
    # plot histogram of R^2 values
    # plt.suptitle(stitles[idx], loc='left')
    ax = plt.subplot(121)
    plt.title(stitles[idx], loc='left')
    format_panels(ax)
    # plt.title(r'$\mathrm{R^2}$ over resamples', y=1)
    tfish, tcad, rfish, rcad = rsq_pair
    cfish, bfish = np.histogram(rfish, bins=nbins)
    cfish = cfish.astype(float) / 1e4
    plt.step(bfish[:-1], cfish, color=colors[0], alpha=0.75, lw=0.75,
            label='phastCons')
    ccadd, bcadd = np.histogram(rcad, bins=nbins)
    ccadd = ccadd.astype(float) / 1e4
    plt.step(bcadd[:-1], ccadd, color=colors[1], alpha=0.75, lw=0.75,
            label='CADD')
    # plt.hist(rfish, bins=nbins, histtype='step', color=colors[0],
    #          alpha=0.75, lw=0.75, label='phastCons', normed=True)
    # plt.hist(rcad, bins=nbins, histtype='step', color=colors[1],
    #          label='CADD', alpha=0.75, lw=0.75, normed=True)
    # plot true values from each
    plt.axvline(x=tfish, color=colors[0], ls='--', lw=1)
    rstr = r'$\mathrm{R^2}=%.3f$' % tfish
    plt.text(tfish, tpos[idx], rstr, color=colors[0], va='bottom', ha='right',
             rotation=90, fontsize=8)
    plt.axvline(x=tcad, color=colors[1], ls='--', lw=1)
    rstr = r'$\mathrm{R^2}=%.3f$' % tcad
    plt.text(tcad+tplus[idx], tpos[idx], rstr, color=colors[1], va='bottom',
             ha='left', rotation=90, fontsize=8)
    # plt.xticks(rotation=90)
    # plt.yticks([100, 200, 300, 400], x=0.04)
    # plt.ylabel(r'# resamples', labelpad=2)
    # plt.ylim(0, ymax)

    # RELATIVE FREQUENCY VERSION
    ymax /= 1e4
    ytcks = np.array([100, 200, 300, 400]) / 1e4
    plt.yticks(ytcks, x=0.04)
    plt.ylabel(r'relative frequency', labelpad=2)
    plt.ylim(0, ymax)
    plt.xticks(y=0.04)

    if idx == 2:
        plt.xticks([0.54, 0.57, 0.6, 0.63], y=0.04)
    # plt.xlim(*xl1[idx])
    plt.xlabel(r'resampled $\mathrm{R^2}$', labelpad=1)
    plt.legend(loc='upper left', handlelength=0.5, labelspacing=0.25,
               borderpad=-0.1)

    # plot histogram of differences in R^2
    ax = plt.subplot(122)
    format_panels(ax)
    # plt.title(r'ii) $\mathrm{R^2_{phastCons}}-\mathrm{R^2_{CADD}}$', y=1)
    cnts, bins = np.histogram(rfish-rcad, bins=nbins)
    cnts = cnts.astype(float) / 1e4
    # cnts, bins = np.histogram(rcad-rfish, bins=nbins)
    widths = bins[1:] - bins[:-1]
    xi = bins[:-1] + widths/2.0
    mneg = (xi<0)
    mpos = (xi>=0)
    plt.bar(xi[mneg], cnts[mneg], widths[mneg], color=colors[1], alpha=0.75)
    plt.bar(xi[mpos], cnts[mpos], widths[mpos], color=colors[0], alpha=0.75)
    plt.axvline(0, color='k', ls=':', lw=0.5)
    # pval = 1.0*np.sum(rfish-rcad>0)/len(rcad)
    # r_plus1 = 1.0*np.sum(rcad-rfish>0)+1
    r_plus1 = 1.0 * np.sum(rfish-rcad > 0) + 1
    n_plus1 = len(rcad) + 1
    pval = r_plus1 / n_plus1

    pstr = r'$p={:.3f}$'.format(pval)
    plt.text(0.855, 0.8, pstr, transform=plt.gcf().transFigure, fontsize=10)

    # plt.yticks([100, 200, 300, 400], x=0.04)
    # plt.ylim(0, ymax)
    # plt.ylabel(r'# resamples', labelpad=2)
    plt.yticks(ytcks, x=0.04)
    plt.ylabel(r'relative frequency', labelpad=2)
    plt.ylim(0, ymax)
    plt.xticks(y=0.04)
    # plt.xticks(xt2[idx], y=0.04)
    # plt.xlim(*xl2[idx])
    # plt.xlabel(r'$\mathrm{R^2}$ difference', labelpad=2)
    # plt.xlabel(r'$\mathrm{R^2_{CADD}}-\mathrm{R^2_{phastCons}}$', labelpad=2)
    plt.xlabel(r'$\mathrm{R^2_{phastCons}}-\mathrm{R^2_{CADD}}$', labelpad=2)

    # plt.legend()
    fsave = final_dir + '/sfigs/updated_Mar2021.bootstrap.rsq.{}.png'.format(wins[idx])
    plt.savefig(fsave, dpi=512)
    plt.close()


for idx in xrange(3):
    plot_rsq_bootstraps_2(bootstraps[idx], idx)


#%% RESAMPLE FOR R^2 BOOTSTRAP (FISH VS APE) -- OLD, NOT USED
def resample_maps_2(win, n_resamples):
    """resample from main best models and calculate R^2"""
    # get both sets of results
    obs1, prd1 = get_pop_array(win, 'YRI', 'cadd94_gmask', 'self')
    obs2, prd2 = get_pop_array(win, 'YRI', 'cadd98_gmask', 'self')

    # set size of random integer space to sample data indices from
    imax = len(obs1)

    # PATCH
    if len(obs1) > len(obs2):
        imax = len(obs2)
        size_diff = len(obs1) - len(obs2)
        msk = np.ones(shape=len(obs1), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs1))
        msk[rnd_rm] = False
        print size_diff, rsquared_function(obs1, prd1)
        obs1, prd1 = obs1[msk], prd1[msk]
        print rsquared_function(obs1, prd1)
    if len(obs2) > len(obs1):
        size_diff = len(obs2) - len(obs1)
        msk = np.ones(shape=len(obs2), dtype=bool)
        rnd_rm = randint_unique(size_diff, len(obs2))
        msk[rnd_rm] = False
        print size_diff, rsquared_function(obs2, prd2)
        obs2, prd2 = obs2[msk], prd2[msk]
        print rsquared_function(obs2, prd2)

    assert len(obs1) == len(obs2)

    # get the true R^2 for each run
    true_1 = rsquared_function(obs1, prd1)
    true_2 = rsquared_function(obs2, prd2)

    # lists for R^2 results from bootstraps
    rsq_1, rsq_2 = [], []

    # # set size of random integer space to sample data indices from
    # imax = len(obs1)

    # do N resamples of each data set anc calculate the R^2 values
    for _ in xrange(n_resamples):
        # sample indices resample with replacement
        idx = np.random.randint(0, imax, size=imax)
        # calculate R^2 from each set of data using the same idx
        r1 = rsquared_function(obs1[idx], prd1[idx])
        rsq_1.append(r1)
        r2 = rsquared_function(obs2[idx], prd2[idx])
        rsq_2.append(r2)

    return true_1, true_2, np.array(rsq_1), np.array(rsq_2)


# tfish, tcad, rfish, rcad = resample_maps(1e6, 10000)
bstp_fish_ape = [resample_maps_2(w, 10000) for w in [1e6]]


#%% BOOTSTRAP PLOTS COMPARING R^2 (FISH/CADD VS APE) -- OLD, NOT USED
def plot_rsq_bootstraps_3(rsq_pair, idx, nbins=100):
    """plot the distribution of boostrap R^2 values and their differences"""
    # idx += 2
    stitles = ['A) 62.5kb windows', 'B) 250kb windows', 'C) 1Mb windows']
    # stitles = ['A) 10kb windows', 'B) 100kb windows', 'C) 1Mb windows']
    # wins = ['10kb', '100kb', '1Mb']
    wins = ['63kb', '250kb', '1Mb']
    colors = ['darkorange', 'firebrick']
    lbls = ['CADD 6%', 'CADD 2%']
    xt1 = [[0.145, 0.15, 0.155, 0.16], [0.38, 0.4, 0.42], [0.55, 0.6, 0.65]]
    xl1 = [(0.1455, 0.1595), (0.365, 0.425), (0.54, 0.66)]
    xt2 = [[-0.008, 0, 0.008], [-0.02, 0, 0.02], [-0.01, 0, 0.01]]
    xl2 = [(-0.011, 0.011), (-0.03, 0.03), (-0.016, 0.016)]
    tpos = [5, 5, 5]
    tplus = [0.0006, 0.0015, 0.003]
    plt.figure(figsize=(4.32, 2.16))
    plt.subplots_adjust(left=0.12, top=0.78, right=0.995, wspace=0.2,
                        bottom=0.17, hspace=0.4)
    # plot histogram of R^2 values
    plt.suptitle(stitles[idx])
    ax = plt.subplot(121)
    format_panels(ax)
    plt.title(r'i) $\mathrm{R^2}$ over resamples', y=1)
    tfish, tcad, rfish, rcad = rsq_pair
    plt.hist(rfish, bins=nbins, histtype='step', color=colors[0],
             alpha=0.75, lw=0.75, label=lbls[0])
    plt.hist(rcad, bins=nbins, histtype='step', color=colors[1],
             label=lbls[1], alpha=0.75, lw=0.75)

    # plot true values from each
    plt.axvline(x=tfish, color=colors[0], ls='--', lw=1)
    rstr = r'$\mathrm{R^2}=%.3f$' % tfish
    txt = plt.text(tfish+tplus[idx], 5, rstr, color=colors[0], va='bottom',
                   ha='left', rotation=90, fontsize=8, weight='bold')
    # bbox=dict(facecolor='white', alpha=0.5, lw=0)
    txt.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='k'),
                           path_effects.Normal()])
    plt.axvline(x=tcad, color=colors[1], ls='--', lw=1)
    rstr = r'$\mathrm{R^2}=%.3f$' % tcad
    txt = plt.text(tcad, tpos[idx], rstr, color=colors[1], va='bottom',
                   ha='right', rotation=90, fontsize=8, weight='bold')
    txt.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='k'),
                           path_effects.Normal()])
    plt.yticks([100, 200, 300, 400], x=0.04)
    plt.ylabel(r'# resamples', labelpad=2)
    plt.ylim(0, 450)
    plt.xticks(y=0.04)
    # plt.xlim(*xl1[idx])
    plt.xlabel(r'resample $\mathrm{R^2}$', labelpad=1)
    plt.legend(loc='upper left', handlelength=0.5, labelspacing=0.25,
               borderpad=-0.1)

    # plot histogram of differences in R^2
    ax = plt.subplot(122)
    format_panels(ax)
    plt.title(r'ii) $\mathrm{R^2_{6\%}}-\mathrm{R^2_{2\%}}$', y=1)
    cnts, bins = np.histogram(rfish-rcad, bins=nbins)
    widths = bins[1:] - bins[:-1]
    xi = bins[:-1] + widths/2.0
    mneg = (xi<0)
    mpos = (xi>=0)
    plt.bar(xi[mneg], cnts[mneg], widths[mneg], color=colors[1], alpha=0.75)
    plt.bar(xi[mpos], cnts[mpos], widths[mpos], color=colors[0], alpha=0.75)
    plt.axvline(0, color='k', ls=':', lw=0.5)
    pval = 1.0*np.sum(rfish-rcad>0)/len(rcad)
    pstr = r'$p={:.3f}$'.format(pval)
    plt.text(0.81, 0.65, pstr, transform=plt.gcf().transFigure)

    plt.yticks([100, 200, 300, 400], x=0.04)
    plt.ylim(0, 450)
    plt.xticks(y=0.04)
    # plt.xlim(*xl2[idx])
    plt.xlabel(r'$\mathrm{R^2}$ difference', labelpad=2)
    # plt.legend()
    fsave = final_dir + '/sfigs/bootstrap.rsq.fish_ape.{}.png'.format(wins[idx])
    plt.savefig(fsave, dpi=512)
    plt.close()


# for idx in xrange(3):
#     plot_rsq_bootstraps_3(bstp_fish_ape[idx], idx)
plot_rsq_bootstraps_3(bstp_fish_ape[0], 2)


#%%