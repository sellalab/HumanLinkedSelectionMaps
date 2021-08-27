__author__ = 'davidmurphy'


import os
import seaborn
import tarfile
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir, \
    human_autosomes, RunStruct, cst_from_fldr
from classes.phylotree import parse_exptotsub
from data_processing.functions import rsquared_function
from data_processing.data_tools import cg_hypermutable
from figures.common_functions import get_bbins


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'
cst = ChromStruct('chr1')
cg_dict = cg_hypermutable()
#%%
bstr = final_dir + '/cadd93_bth600/bmap/{}.bmap.txt'
nstr = root_dir + '/data/nmsk/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'
cg_bvals = []
cg_neut = []
all_bv, all_lv = [], []
all_neut = []
for ch in human_autosomes:
    cst.chrom = ch
    nmsk = np.load(nstr.format(ch))['neutmask']
    barr = np.zeros(shape=cst.chlen)
    bv, lv = np.loadtxt(bstr.format(ch)).T
    start = np.concatenate(([0], np.cumsum(lv[:-1]))).astype(int)
    for i in xrange(len(start)):
        all_neut.append(np.sum(nmsk[start[i]:start[i]+int(lv[i])]))
        barr[start[i]:start[i]+int(lv[i])] = bv[i]
        all_bv.append(bv[i])
        all_lv.append(lv[i])
    for (si, sj) in cg_dict[ch]:
        cg_neut.append(np.sum(nmsk[si:sj]))
        bmean = np.mean(barr[si:sj][nmsk[si:sj]])
        cg_bvals.append(bmean)
#%%
f_cg = root_dir + '/data/coords/CG_hypermutable.txt'
rel_pos = []
with open(f_cg, 'r') as f:
    for line in f:
        ch, pos, stat = line.strip('\n').split()
        if stat == '1':
            cst.chrom = ch
            rpos = float(pos) * 1e6 / cst.chlen
            rel_pos.append(rpos)

#%%
f_save = final_dir + '/sfigs/cg_hyper_spatial_dist.png'
flat_weights = (1.0 / len(rel_pos)) * np.ones(shape=len(rel_pos))
plt.figure(figsize=(6.5, 3.25))
plt.subplots_adjust(bottom=0.15, left=0.1, top=1, right=1)
plt.hist(rel_pos, weights=flat_weights, bins=25)
plt.ylabel('fraction of CG hypermutable tracts')
plt.xlabel('relative chromosomal position')
plt.savefig(f_save, dpi=512)
plt.close()
#%%
cgmsk = ~np.isnan(cg_bvals)
cg_bvals = np.array(cg_bvals)
cg_neut = np.array(cg_neut)
f_save = final_dir + '/sfigs/cg_hyper_bval_dist.png'
plt.figure(figsize=(6.5, 3.25))
plt.subplots_adjust(bottom=0.15, left=0.1, top=1, right=1)
plt.hist(cg_bvals[cgmsk], weights=cg_neut[cgmsk], bins=50, histtype='step', label='CG hypermutable B vals',
         normed=1, lw=1.25)
plt.hist(all_bv, weights=all_neut, bins=50, histtype='step', label='All B vals',
         normed=1, lw=1.25)
plt.ylabel('fraction of sites')
plt.xlabel('B value')
plt.legend(loc='upper left')
plt.savefig(f_save, dpi=512)
plt.close()
#%%
tkn = 'cadd93_extel.filter.gmap.edge.0.1'
bbin = get_bbins(tkn)
#%%
nstr = root_dir + '/data/nmsk/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'
cg_arrs = {}
cst = ChromStruct('chr1')
for ch in human_autosomes:
    cst.chrom = ch
    nmsk = np.load(nstr.format(ch))['neutmask']
    c_arr = np.zeros(shape=cst.chlen)
    for (si, sj) in cg_dict[ch]:
        c_arr[si:sj] = 1
    # turn non-neutral sites to 0s
    c_arr[~nmsk] = 0
    cg_arrs[ch] = c_arr
#%%
cg_counts = []
for (i, b) in enumerate(bbin):
    n_sites = 0
    for (c, start, end, num) in b:
        ch = 'chr{}'.format(int(c))
        start, end = int(start), int(end)
        n_cg = cg_arrs[ch][start:end].sum()
        if n_cg > num:
            print i, n_cg, ch, start, end, num
            break
        n_sites += n_cg
    cg_counts.append(n_sites)
#%%
pi = np.arange(100)
flat = np.ones(shape=100)*0.01
cg = np.array(cg_counts)
f_save = final_dir + '/sfigs/cg_hyper_b_percentile.png'
plt.figure(figsize=(6.5, 3.25))
plt.subplots_adjust(bottom=0.15, left=0.12, top=1, right=1)
plt.axhline(0.01, label='all neutral', color='k')
plt.scatter(pi, cg / np.sum(cg, dtype=float), color='r', label='CG hypermutable')
plt.ylabel('fraction of neutral sites per bin')
plt.xlabel('B percentile')
plt.ylim(0, 0.022)
plt.legend(loc='lower right')
plt.savefig(f_save, dpi=512)
plt.close()
#%%
fdir = root_dir + '/result/final_files/{}/'.format('cadd93_bth600')
f_sort = fdir + 'sort_gc_cm_cn_n{}.txt'.format(100)
pred = np.loadtxt(f_sort)[:,2]
pred /= pred.max()
f_save = final_dir + '/sfigs/cg_hyper_pred_percentile.png'
plt.figure(figsize=(6.5, 3.25))
plt.subplots_adjust(bottom=0.15, left=0.12, top=1, right=1)
plt.axhline(0.01, label='all neutral', color='k')
plt.scatter(pred, cg / np.sum(cg, dtype=float), color='r', label='CG hypermutable')
plt.ylabel('fraction of neutral sites per bin')
plt.xlabel('prediction')
plt.ylim(0, 0.022)
plt.legend(loc='lower right')
plt.savefig(f_save, dpi=512)
plt.close()
#%%