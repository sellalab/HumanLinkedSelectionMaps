#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 23:29:52 2018

@author: davidmurphy
"""

import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as pcor
from precalc.lh_inputs import load_saved
from classes.runstruct import ChromStruct, root_dir, izip
from data_processing.functions import rsquared_function, swap_root
from likelihood.cllh_functions import predicted_pi


#%%
def plot_sorted(sarr, cst, num_bins):
    # set plot title
    title = 'prediction sorted pi'

    # expand sorted array results
    div, pi, pr = sarr.T

    # TODO: get pi0 properly
    meanpi0 = pr.max()
    meanpi = cst.stat.meanpi / meanpi0
    x = np.arange(len(div))
    pi /= meanpi0
    pr /= meanpi0

    corr, pval = pcor(pi, pr)

    gray = 'DarkSlateGray'
    fig = plt.figure(figsize=(10, 7))
    plt.subplots_adjust(left=0.13, bottom=0.12, right=0.78, top=0.90, wspace=0.05,
                        hspace=0.2)
    fig.add_subplot(111)

    plt.plot(x, pi, label='observed diversity', color=gray, lw=3)
    plt.plot(x, pr, label='predicted', color='Fuchsia', alpha=0.75, lw=3)

    message = ' {}\n {:>13}'.format(*'no background;selection'.split(';'))
    plt.axhline(y=1, color='k', ls='--', lw=2.0)
    # plt.text((0.5 * num_bins), 1.005, 'w/o background selection', ha='center', va='bottom', fontsize=18)
    plt.text((1.025 * num_bins), 1, message, ha='left', va='center', fontsize=20)

    message = '   {:>6}\n   {}'.format(*'mean;diversity'.split(';'))
    plt.axhline(y=meanpi, color=gray, ls=':', lw=2.0)
    # plt.text((0.8 * num_bins), meanpi - 0.005, 'mean diversity', ha='center', va='top', fontsize=18, color=gray)
    plt.text((1.025 * num_bins), meanpi, message, ha='left', va='center',
             fontsize=20, color=gray)

    # Pearson correlation
    plt.text((0.5 * num_bins), 0.55, r'$Pearson\ R^2 = {:.4f}$'.format(corr),
             ha='center', va='center', fontsize=20)

    # title w/ species or cons info
    plt.title(title, fontsize=20)
    # plt.title(rst.token, fontsize=20)
    # plt.title(init_file.split('/')[-1].split('.')[0], fontsize=20)

    plt.xlabel('strength of background selection', fontsize=20)
    plt.xticks(fontsize=18)
    plt.xlim((-0.01 * num_bins), (1.01 * num_bins))

    plt.ylabel('scaled diversity', fontsize=20, labelpad=10)
    plt.yticks(fontsize=18)
    plt.ylim(0.37, 1.23)

    plt.legend(prop={'size': 18}, loc='upper left', ncol=1, numpoints=3,
               borderaxespad=0.1)
    # fsave = '/Users/davidmurphy/Desktop/sort.fig.png'
    fsave = '/Users/davidmurphy/Desktop/unbound.sort.fig.png'

    plt.savefig(fsave, dpi=256)
    # plt.show()


def calc_rsquared(scale, cum_pos, ms, nt, dv, pred):
    # create genomic windows for given scale
    windows = np.arange(0, cum_pos[-1], scale)
    # create upper, lower indices to sort data into genomic windows
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    obs, prd, num = [], [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate average pi for the window
            pi = nt[i:j, 1].sum() / nt[i:j].sum()
            # calculate average predicted pi for window (weighted by sites)
            pr = np.average(pred[i:j], weights=dv[i:j])
            # calculate site count for window
            n = dv[i:j].sum()
            # record stats to lists
            obs.append(pi)
            prd.append(pr)
            num.append(n)

    # convert all results to arrays
    obs = np.array(obs)
    prd = np.array(prd)
    num = np.array(num)
    msk = (num > 0.2*np.mean(num))  # TODO: what to do about this filter??

    # calculate R squared
    rsq = rsquared_function(obs[msk], prd[msk])

    return rsq


def basic_sort(neut, div, pred, num):
    """
    Sort data arrays by predictions and average results into bins.
    :type init: RunStruct innit file
    :param num: the number of bins to sort into
    """
    # get sorting indices based on prediction value
    sidx = np.argsort(pred)
    # sort predictions and div/poly data
    neut, div, pred = [a[sidx].astype('f8') for a in neut, div, pred]
    sites, subs = div.T
    # calculate mean div
    meandiv = subs.sum() / sites.sum()

    # get start and end indices per partition
    idx = sortbin_edges(sites=sites, numbins=num)

    # gather data summaries
    sorted_array = []
    for (i, j) in idx:
        # calculate each statistic for the present bin
        norm_div = np.sum(subs[i:j]) / (meandiv * np.sum(sites[i:j]))
        scaled_pi = np.sum(neut[i:j, 1]) / (np.sum(neut[i:j]) * norm_div)
        # use weighted mean for pred
        prediction = np.average(pred[i:j], weights=sites[i:j])
        sorted_array.append([norm_div, scaled_pi, prediction])

    return np.array(sorted_array)


def sortbin_edges(sites, numbins):
    """
    get the upper indices of sorted data array that divide
    data into bins of equal neutral site counts
    """

    # get number of sites needed so that numbins x numsites = total sites
    numsites = int(np.sum(sites) / numbins)

    # find indices that partition cumulative sorted site count into numsites
    cumsites = np.cumsum(sites)
    bounds = np.arange(numsites, cumsites[-1], numsites)

    # get the ending index for each partition
    jdx = list(np.searchsorted(a=cumsites, v=bounds))

    # return a list of (start, end) indices for each partition
    return zip([0] + jdx[:-1], jdx)


#%%
# directory to result files
# folder_name = 'std_run_feb2019'
folder_name = 'unboundB_run_feb2019'
fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
# list of result files
flst = [fdir+f for f in os.listdir(fdir) if f.endswith('.txt')]
# list of chromstructs for each result
rlst = [ChromStruct('chr1', init=f) for f in flst]

# create tuples of file index, final LH result
ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
# get indices of the top 3 best LH runs (after sorting on LH)
top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]
# take the average of the params of the top 3 runs
pavg = np.average([rlst[i].params for i in top3], axis=0)

# use any file for init
init = flst[0]
swap_root(init)
cst = ChromStruct(chrom='chr1', init=init)

#%%
# create the plot
plt.figure(figsize=(10, 5))
plt.subplots_adjust(left=0.15, wspace=0.65, bottom=0.15)
plt.subplot(1, 7, (1, 6))
plt.title('selection parameters', fontsize=18)

# plot results from the top3 indices of runs (averaged)
x = np.arange(6)
s = 0
w = 0.8
pmfs = []
for i in top3:
    r = rlst[i]
    pmfs.append(r.uvec[0])
pmfs = np.array(pmfs)
mpmf = np.average(pmfs, axis=0) * 1e8

plt.bar(x + s, mpmf, w, color='darkorange')
plt.ylabel(r'$\mu_{del} \cdot 10^{-8}$', size=18)
plt.yticks(size=18)
plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)], size=18)
plt.xlabel('deleterious fitness effect', size=18)
plt.legend()

# plot udel
sp3 = plt.subplot(177)
sp3.yaxis.tick_right()
sp3.yaxis.set_label_position('right')
plt.yticks(size=18)
plt.title(r'$\mu_{del} \cdot 10^{-8}$', size=18)
plt.bar(0, sum(mpmf), 0.65, color='darkorange')
plt.xlim(-0.5,0.5)
plt.xticks([])
# plt.show()
# fsave = '/Users/davidmurphy/Desktop/params.fig.png'
fsave = '/Users/davidmurphy/Desktop/unbound.params.fig.png'
plt.savefig(fsave, dpi=256)

#%%
# load all arrays
sg, bs, cs, nu, nt, dv, pl = load_saved(cst)
# get the number of sites per segmet from div array
ns = dv[:,0]
# create a mask for segments with no data, filter div array right away
ms = (ns > 0)
# mask additional arrays
nt, ns, nu, bs = [ar[ms] for ar in nt, ns, nu, bs]
# convert bs to 0-1 scale
bs *= np.log1p(-1.0 / cst.bscl)
# apply min filtering to bs values (if optimization used filter)
if folder_name == 'std_run_feb2019':
    bs = np.maximum(np.log(0.01), bs)

#%%
# create constant u array
nu_const = np.full(shape=nu.size, fill_value=cst.stat.meandiv)
# create predictions with constant u
pred_const = predicted_pi(pavg, cst, nu_const, bs, None)
# set the number of bins to average sorted data into
num_bins = 100
# get prediction sorted mean values
sarr = basic_sort(nt, dv, pred_const, num_bins)
plot_sorted(sarr, cst, num_bins)

#%%
# convert segments into positions
cum_pos = np.cumsum(sg)
# calculate predicted pi from arrays and average of top 3 params
pred = predicted_pi(pavg, cst, nu, bs, None)

#%%
fmcv = '{}/result/rsq_maps/old_data_rsq/phylo-2/mcvicker_map_BS1_CS0_161129211053_rsqmap.txt'.format(root_dir)
mc = np.loadtxt(fmcv)

#%%
scales = mc[:,0]
rsq = []
prm = rlst[top3[0]].params
for sc in scales:
    pred = predicted_pi(prm, cst, nu, bs, None)
    rsq.append(calc_rsquared(sc, cum_pos, ms, nt, ns, pred))
    print rsq[-1]

#%%
plt.figure(figsize=(8,8))
plt.subplots_adjust(left=0.16, bottom=0.15, right=0.95, top=0.95, wspace=0,
                    hspace=0)

markprop = dict(markersize=10, markeredgewidth=0, markeredgecolor='1',
                alpha=0.9, linewidth=2, linestyle='None')

xvals = np.log10(mc[:,0])
plt.plot(xvals, mc[:,1], marker='s', label='McVicker prediction',
         color='darkturquoise', **markprop)
plt.plot(xvals, rsq, marker='o', label='Our prediction',
         color='fuchsia', **markprop)
# x-labels
plt.xlabel('window size', fontsize=24, labelpad=10)
xticks = np.arange(4, 6.1, 0.5)
plt.xticks(xticks, ['10kb', '', '100kb', '', '1Mb'], fontsize=22)
plt.xlim(4, 6.35)
# y-labels
plt.ylabel('variance explained', fontsize=24, labelpad=10)
plt.yticks(np.arange(0.1, 0.6, 0.05), fontsize=22)
plt.ylim(0.0675, 0.59)
legprop = dict(size=18)

plt.legend(loc='upper left', prop=legprop, numpoints=1)
# fsave = '/Users/davidmurphy/Desktop/rsq.fig.png'
fsave = '/Users/davidmurphy/Desktop/unbound.rsq.fig.png'
plt.savefig(fsave, dpi=256)
# plt.show()


#%%
scale = 1e6
# create genomic windows for given scale
windows = np.arange(0, cum_pos[-1], scale)
# create upper, lower indices to sort data into genomic windows
jdx = np.searchsorted(cum_pos[ms], windows)
idx = np.concatenate(([0], jdx[:-1]))

cst.chrom = 'chr1'
msk = cum_pos[ms][idx] < cst.chlen

# prepare empty lists for each relevant statisic
pos, obs, prd = [], [], []
for (i, j) in izip(idx[msk], jdx[msk]):
    if j > i:
        # calculate average pi for the window
        pi = nt[i:j, 1].sum() / nt[i:j].sum()
        # calculate average predicted pi for window (weighted by sites)
        pr = np.average(pred[i:j], weights=ns[i:j])
        # record stats to lists
        obs.append(pi)
        prd.append(pr)

    else:
        obs.append(np.nan)
        prd.append(np.nan)

    pos.append(cum_pos[ms][j])


#%%
obs, prd = np.array(obs), np.array(prd)
plt.figure(figsize=(13, 4.25))
plt.subplots_adjust(left=0.2, bottom=0.18, right=0.95, top=0.9, wspace=0.05,
                    hspace=0.1)
plt.plot(obs, color='DarkSlateGray', label='Observed', lw=2.5)
plt.plot(prd, color='Fuchsia', alpha=0.8, label='Our map', lw=2.5)

# plt.title('chromosome 1', fontsize=28)
plt.xlabel('chr1 position (Mb)', fontsize=24)
plt.xticks(range(25, 226, 50), range(25, 226, 50), fontsize=22)
plt.xlim(0, 250)
plt.ylabel('scaled diversity', fontsize=24)
plt.yticks(fontsize=22)
plt.legend(prop={'size': 22}, ncol=3, loc='upper center')
# fsave = '/Users/davidmurphy/Desktop/chr1.fig.png'
fsave = '/Users/davidmurphy/Desktop/unbound.chr1.fig.png'
plt.savefig(fsave, dpi=256)
plt.show()

