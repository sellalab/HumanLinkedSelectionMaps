__author__ = 'davidmurphy'

import os
# import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir

final_dir = root_dir + '/result/final_files'


#%%
def get_bmap(bfile):
    """load bmap and return as pos, b"""
    b, l = np.loadtxt(bfile).T
    b = np.exp(b * np.log1p(-1.0 / 100.0))
    return b, np.cumsum(l)


def get_cmap(cfile):
    """load cmap and return as pos, c"""
    c, l = np.load(cfile)['cvals'].T
    return c, np.cumsum(l)


def get_pred(pfile):
    """load results for chr1 for ape_cons94"""
    prd, obs, num = np.loadtxt(pfile).T
    return prd, obs, num


#%%
# bf = 'AA_Map_10kb.chr1.ape_cons94_clean_extel.t0.00100000.merge1e07.bkgd'
bf = 'AA_Map_10kb.chr1.cds.t0.00010000.bkgd'
bfile = final_dir + '/mainfigs/precalc_maps/' + bf
b, xb = get_bmap(bfile)

cf = 'AA_Map_10kb.chr1.YRI_nonsyn_s1.s0.00031623.npz'
cfile = final_dir + '/mainfigs/precalc_maps/' + cf
c, xc = get_cmap(cfile)
cdenom = 1.0 + c
cv = 1.0 / cdenom
#%%
bi = np.interp(xc, xb, b)**(1.0/2)
bdenom = 1.0 / bi
bv = 1.0 / bdenom
rdenom = bdenom+c
rv = 1.0 / rdenom


#
plt.figure(figsize=(4,2))
plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
xlo, xhi = 1.15e7, 1.2e7
# xlo, xhi = 5.37e7, 5.42e7
msk = (xc>=xlo) & (xc<=xhi)
ax1 = plt.subplot(111)
gridargs = dict(color='k', linestyle=':', lw=0.2, alpha=0.5)
ax1.set_facecolor('white')
ax1.grid(**gridargs)

plt.plot(xc[msk], bv[msk], lw=1.25, color='red')
plt.plot(xc[msk], cv[msk], lw=1.25, color='blue')
plt.plot(xc[msk], rv[msk], lw=1.25, color='purple')

idx = range(0, msk.sum(), 10)
sv = rv[msk][idx] + np.random.normal(0, 0.05, size=len(idx))
plt.plot(xc[msk][idx], sv, lw=1.25, color='gray', alpha=0.5)
# plt.xticks([])
# plt.yticks([])
plt.show()


#%% EXAMPLE BS/CS MAP FRAGMENTS
# plt.figure(figsize=(50, 2))
# xlo, xhi = 5.37e7, 5.42e7
# xlo, xhi = 1.e7, 1.1e7

plt.figure(figsize=(2.5, 1.25))

plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
# msk = (xc>=xlo) & (xc<=xhi)
ax1 = plt.subplot(111)
ax1.set_facecolor('white')
ax1.axis('on')
ax1.patch.set_edgecolor('black')
plt.plot(xc[msk], bv[msk], lw=1.25, color='red')
fsave = final_dir + '/mainfigs/precalc_maps/chr1_bs.png'
plt.xticks([])
plt.yticks([])
plt.ylim(0, 1.1)
plt.savefig(fsave, dpi=512)
plt.close()
#
plt.figure(figsize=(2.5, 1.25))
plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
ax1 = plt.subplot(111)
ax1.set_facecolor('white')
ax1.axis('on')
ax1.patch.set_edgecolor('black')
# xlo, xhi = 6.7e7, 6.8e7
# msk = (xc>=xlo) & (xc<=xhi)
plt.plot(xc[msk], cv[msk], lw=1.25, color='blue')
plt.xticks([])
plt.yticks([])
plt.ylim(0, 1.1)
fsave = final_dir + '/mainfigs/precalc_maps/chr1_cs.png'
plt.savefig(fsave, dpi=512)
plt.close()
#%%
plt.figure(figsize=(2.5, 1.25))
plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
ax1 = plt.subplot(111)
ax1.set_facecolor('white')
ax1.axis('on')
ax1.patch.set_edgecolor('black')
plt.plot(xc[msk][idx], sv, lw=2, color='darkslategray', alpha=1, label='observed')
plt.plot(xc[msk], rv[msk], lw=2, color='darkorange', alpha=0.8, label='predicted')

plt.xticks([])
plt.yticks([])
plt.ylim(-0.2,1.1)
plt.legend(loc='lower center', ncol=2)
fsave = final_dir + '/mainfigs/precalc_maps/chr1_bs+cs.png'
plt.savefig(fsave, dpi=512)
# plt.close()


#%% EXAMPLE PREDICTED/OBSERVED MAP FRAGMENT
pf = 'chr1.1.00e+06win_0.5slide_pred_and_obs.txt'
pfile = final_dir + '/ape_cons94_bth400/' + pf
prd, obs, num = get_pred(pfile)
xi = np.arange(0, prd.size / 2.0, 0.5)
pi_mean = np.nanmean(obs)
prd /= pi_mean
obs /= pi_mean
obs[(obs < 0.25) | (obs > 1.75)] = np.nan


#%%
xlo, xhi = 160, 200
msk = (xi>=xlo) & (xi<=xhi)
plt.figure(figsize=(4, 2))
plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
plt.plot(xi[msk], obs[msk], color='gray', lw=2, label='observed')
plt.plot(xi[msk], prd[msk], color='purple', lw=2, label='predicted')
plt.xticks([])
plt.yticks([])
plt.legend()
# plt.ylim(0, 1.1)
fsave = final_dir + '/mainfigs/precalc_maps/chr1_pred.png'
plt.savefig(fsave, dpi=512)
plt.close()


#%% EXAMPLE DISTRIBUTION OF FITNESS EFFECTS
udel = [0.2, 0.1, 0.7, 1.2, 0.1, 0.3]
alph = [0.1, 0.5, 0.2, 0.1, 1, 0.2]
xbar = np.arange(6)
plt.figure(figsize=(2.5, 1.25))
plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
plt.bar(xbar, udel, color='red')
plt.xticks([])
plt.yticks([])
# plt.ylim(0, 1.1)
fsave = final_dir + '/mainfigs/precalc_maps/bs_params.png'
plt.savefig(fsave, dpi=512)
plt.close()

plt.figure(figsize=(2.5, 1.25))
plt.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
plt.bar(xbar, alph, color='blue')
plt.xticks([])
plt.yticks([])
# plt.ylim(0, 1.1)
fsave = final_dir + '/mainfigs/precalc_maps/cs_params.png'
plt.savefig(fsave, dpi=512)
plt.close()


#%%