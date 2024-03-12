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
from skmisc import loess


# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'

#%%

def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def get_sub_rates_2(fldr):
    rst = cst_from_fldr(fldr)
    bbin_dir = root_dir + '/data/phast/bbins'
    res = []
    msk = []
    for bbin in range(1500, 2000):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        # f_name = '{}/bbin{}.exptotsub'.format(bbin_dir, bbin)
        f_name = '{}/{}.bbin{}.exptotsub'.format(bbin_dir, rst.tkn, bbin)
        if not os.path.isfile(f_name):
            print '{} missing bbin {}'.format(fldr, bbin)
            msk.append(False)
            res.append([1]*8)

            continue

        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m
        # adjust total for number of branches
        for i in range(4):
            mat[i,i] /= 14.0
        # get A>C/T>G
        at_tot = (mat[0, :].sum() + mat[3, :].sum())
        # at_tot = (mat[0, 0].sum() + mat[3, 3].sum())

        actg = (mat[0, 1] + mat[3, 2]) / at_tot
        # get A>G/T>C
        agtc = (mat[0, 2] + mat[3, 1]) / at_tot
        # get A>T/T>A
        atta = (mat[0, 3] + mat[3, 0]) / at_tot
        # get C>A/G>T
        cg_tot = (mat[1, :].sum() + mat[2, :].sum())
        # cg_tot = (mat[1, 1].sum() + mat[2, 2].sum())
        cagt = (mat[1, 0] + mat[2, 3]) / cg_tot
        # get C>G/G>C
        cggc = (mat[1, 2] + mat[2, 1]) / cg_tot
        # get C>T/G>A
        ctga = (mat[1, 3] + mat[2, 0]) / cg_tot

        row = actg, agtc, atta, cagt, cggc, ctga, at_tot, cg_tot
        res.append(row)
        msk.append(True)

    res = np.array(res)

    dv_vals = []
    for r in res[:,:6].T:
        dv = r/r[msk].mean()
        dv_vals.append(dv)
    dv = np.array(dv_vals)

    return res, dv, np.array(msk)


# LOAD *FIXED* SNP DATA
def get_snp_poly(fldr, nbins=100):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    # f_sort = fdir + 'predsort_more.txt'
    # f_sort = fdir + 'basic_sort_n{}.txt'.format(nbins)
    f_sort = fdir + 'sort_gc_cm_cn_n{}.txt'.format(nbins)

    f_ssnp = fdir + 'sort_snptype_fix_anc_n{}.txt'.format(nbins)

    snp = np.loadtxt(f_ssnp)
    # div, pi, pred, gc, cn, cm = np.loadtxt(f_sort).T
    div, pi, pred = np.loadtxt(f_sort)[:,:3].T

    if fldr == 'cadd93_extel_rm_CG':
        rst = cst_from_fldr('cadd93_bth600')
    else:
        rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    all_pimean = pi.mean()
    pred_mean = pred.mean()

    # pred *= (all_pimean/(pred_mean*pi0))
    pred /= pi0
    pi_vals = []
    for s in snp[:,:6].T:
        pi = s*all_pimean/(s.mean()*pi0)
        pi_vals.append(pi)
    ancpi = np.array(pi_vals)

    return ancpi, pred, snp


#
msk = get_sub_rates_2('cadd93_extel_rm_CG')[2]
# fl = 'cadd93_bth600'
fl = 'cadd93_extel_rm_CG'
res, dv, _ = get_sub_rates_2(fl)
ancpi, pred, snp = get_snp_poly(fl, nbins=2000)

res = res[msk]
dv = dv.T[msk].T
ancpi = ancpi.T[1500:][msk].T
pred = pred[1500:][msk]
snp = snp[1500:][msk]
print fl, pred.mean()

ymin = 0.65
ymax = 1.55
lspan = 0.2

# AT ORIGINATING 2
use_idx = [0, 1, 2]
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
f_save = root_dir + '/result/final_files/sfigs/{}_AT_2000.png'.format(fl)
plt.figure(figsize=(3, 21.0 / 4.0))
plt.subplots_adjust(top=0.95, right=0.975, hspace=0.1, bottom=0.1, left=0.18)
plt.subplot(111)
if 'CG' in fl:
    plt.title('AT originating (C>G filtered)')
else:
    plt.title('AT originating')
lo_fits = []
for i in use_idx:
    yi = ancpi[i] / dv[i]
    ploess = predict_loess(pred, yi, None, lspan, pred)
    lo_fits.append(ploess)
    plt.plot(pred, yi, marker='o', ms=2, lw=0, alpha=0.8,
             label=labs[i], color=clist[i])

# plot LOESS fits on top (darkened)
for (i, ploess) in zip(use_idx, lo_fits):
    plt.plot(pred, ploess, color=clist[i], lw=2, alpha=0.9)
    # plt.plot(pred, ploess, color='white', alpha=0.5)

plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
# xtick = np.arange(0.5, 1.01, 0.1)
# ytick = np.arange(0.5, 1.41, 0.1)
# plt.xticks(color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
# plt.yticks(ytick)
plt.ylim(ymin, ymax)
plt.xlim(0.9, 0.995)
plt.legend(loc='upper left', ncol=2, frameon=1, framealpha=0.75,
           facecolor='white')
# plt.text(0.01, 0.98, 'A', transform=plt.gcf().transFigure)
#
# plt.subplot(212)
# for i in use_idx:
#     plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
#              label=labs[i], color=clist[i])
# # xtick = np.arange(0.5, 1.01, 0.1)
# # plt.xticks(xtick)
# plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
# plt.xlim(0.9, 1.0)
# plt.ylabel(r'$D/\bar{D}$', labelpad=3)
# # ytick = np.arange(0.9, 1.31, 0.1)
# # plt.yticks(ytick)
# # plt.ylim(0.82, 1.38)
# plt.legend(loc='upper left', frameon=1, framealpha=0.75, facecolor='white')

plt.savefig(f_save, dpi=512)
plt.close()

# GC ORIGINATING 2
use_idx = [3, 4, 5]
clist = ['darkorange', 'red', 'purple', 'fuchsia', 'steelblue', 'rosybrown']
labs = ['A>C/T>G', 'A>G/T>C', 'A>T/T>A', 'C>A/G>T', 'C>G/G>C', 'C>T/G>A']
f_save = root_dir + '/result/final_files/sfigs/{}_GC_2000.png'.format(fl)
plt.figure(figsize=(3, 21.0 / 4.0))
plt.subplots_adjust(top=0.95, right=0.975, hspace=0.1, bottom=0.1, left=0.18)
plt.subplot(111)

if 'CG' in fl:
    plt.title('CG originating (C>G filtered)')
else:
    plt.title('CG originating')

# plot results
lo_fits = []
for i in use_idx:
    yi = ancpi[i] / dv[i]
    ploess = predict_loess(pred, yi, None, lspan, pred)
    lo_fits.append(ploess)
    plt.plot(pred, yi, marker='o', ms=2, lw=0, alpha=0.8,
             label=labs[i], color=clist[i])

# plot LOESS fits on top (darkened)
for (i, ploess) in zip(use_idx, lo_fits):
    plt.plot(pred, ploess, color=clist[i], lw=2, alpha=0.9)
    # plt.plot(pred, ploess, color='white', alpha=0.5)
# plot y=x
plt.plot([0, 1.1], [0, 1.1], label=r'$y=x$', color='darkslategray',
         ls='--', alpha=0.65)
# xtick = np.arange(0.5, 1.01, 0.1)
# ytick = np.arange(0.5, 1.41, 0.1)
# plt.xticks(color='white')
plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)
# plt.yticks(ytick)
plt.ylim(ymin, ymax)
plt.xlim(0.9, 0.995)
plt.legend(loc='upper left', ncol=2, frameon=1, framealpha=0.75,
           facecolor='white')
# plt.text(0.01, 0.98, 'A', transform=plt.gcf().transFigure)

# plt.subplot(212)
# for i in use_idx:
#     plt.plot(pred, dv[i], marker='o', ms=5, lw=0, alpha=0.75,
#              label=labs[i], color=clist[i])
# # xtick = np.arange(0.5, 1.01, 0.1)
# # plt.xticks(xtick)
# plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
# plt.xlim(0.9, 1.0)
# plt.ylabel(r'$D/\bar{D}$', labelpad=3)
# # ytick = np.arange(0.9, 1.31, 0.1)
# # plt.yticks(ytick)
# # plt.ylim(0.82, 1.38)
# plt.legend(loc='upper left', frameon=1, framealpha=0.75, facecolor='white')
# # plt.text(0.01, 0.47, 'B', transform=plt.gcf().transFigure)
plt.savefig(f_save, dpi=512)
plt.close()
#%%