__author__ = 'davidmurphy'

import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir
from figures.other_code.summary_slide import cst_from_fldr

final_dir = root_dir + '/result/final_files'

#%%

def rsq_three_scales(fldr):
    r_file = final_dir + '/{}/rsq.log'.format(fldr)
    w, r = np.loadtxt(r_file).T
    ikeep = np.array([0, 7, 13])

    return r[ikeep]


#%%
ape_rsq = []
cadd_rsq = []

ape_utot = []
cadd_utot = []

#%% COMPARE R^2 AT 1MB
gm_range = range(1, 21)


apefmt = 'ape94_gmask_{:02}'
caddfmt = 'cadd93_gmask_{:02}'

for gm in gm_range:
    # make folders for current gdist
    ape_fldr = apefmt.format(gm)
    cadd_fldr = caddfmt.format(gm)

    # get R^2 at 3 spatial scales
    ape_rsq.append(rsq_three_scales(ape_fldr))
    cadd_rsq.append(rsq_three_scales(cadd_fldr))

    # get utot for ape/cadd
    ape_utot.append(cst_from_fldr(ape_fldr).stat.utot[0])
    cadd_utot.append(cst_from_fldr(cadd_fldr).stat.utot[0])

#%%
ape_rsq = np.array(ape_rsq)
cadd_rsq = np.array(cadd_rsq)
ape_utot = np.array(ape_utot)
cadd_utot = np.array(cadd_utot)

#%%
gi = np.concatenate((np.arange(0, 0.1, 0.01), np.arange(0.1, 2.1, 0.1)))
plt.figure(figsize=(7.5, 2.5))
plt.subplots_adjust(bottom=0.15, left=0.1, right=1, top=1)
plt.plot(gi, ape_rsq[:,0], marker='s', label='ape 10kb', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, ape_rsq[:,1], marker='D', label='ape 100kb', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, ape_rsq[:,2], marker='o', label='ape 1Mb', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_rsq[:,0], marker='s', label='CADD 10kb', color='dodgerblue',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_rsq[:,1], marker='D', label='CADD 100kb', color='dodgerblue',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_rsq[:,2], marker='o', label='CADD 1Mb', color='dodgerblue',
         alpha=0.7, ms=5)
plt.legend()
plt.xlabel('cM of neutral sites removed from telomeres')
plt.ylabel(r'variance explained ($R^2$)')
f_save = final_dir + '/sfigs/rsq_compare_gmasks.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%% R^2 AT 1MB
plt.figure(figsize=(7.5, 2.5))
plt.subplots_adjust(bottom=0.15, left=0.1, right=1, top=1)
plt.plot(gi, ape_rsq[:,2], marker='o', label='ape 1Mb', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_rsq[:,2], marker='o', label='CADD 1Mb', color='dodgerblue',
         alpha=0.7, ms=5)
plt.legend()
plt.xlabel('cM of neutral sites removed from telomeres')
plt.ylabel(r'variance explained ($R^2$)')
f_save = final_dir + '/sfigs/rsq_1Mb_compare_gmasks.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%% R^2 AT 100KB
plt.figure(figsize=(7.5, 2.5))
plt.subplots_adjust(bottom=0.15, left=0.1, right=1, top=1)
plt.plot(gi, ape_rsq[:,1], marker='o', label='ape 100Kb', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_rsq[:,1], marker='o', label='CADD 100Kb', color='dodgerblue',
         alpha=0.7, ms=5)
plt.legend()
plt.xlabel('cM of neutral sites removed from telomeres')
plt.ylabel(r'variance explained ($R^2$)')
f_save = final_dir + '/sfigs/rsq_100kb_compare_gmasks.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%% R^2 AT 10KB
plt.figure(figsize=(7.5, 2.5))
plt.subplots_adjust(bottom=0.15, left=0.1, right=1, top=1)
plt.plot(gi, ape_rsq[:,0], marker='o', label='ape 10Kb', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_rsq[:,0], marker='o', label='CADD 10Kb', color='dodgerblue',
         alpha=0.7, ms=5)
plt.legend()
plt.xlabel('cM of neutral sites removed from telomeres')
plt.ylabel(r'variance explained ($R^2$)')
f_save = final_dir + '/sfigs/rsq_10kb_compare_gmasks.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%% UTOT
plt.figure(figsize=(7.5, 2.5))
plt.subplots_adjust(bottom=0.15, left=0.1, right=1, top=1)
plt.plot(gi, ape_utot/1.4e-8, marker='o', label='ape conserved', color='darkorange',
         alpha=0.7, ms=5)
plt.plot(gi, cadd_utot/1.4e-8, marker='o', label='CADD', color='dodgerblue',
         alpha=0.7, ms=5)
plt.legend()
plt.xlabel('cM of neutral sites removed from telomeres')
plt.ylabel(r'$u_{del}\ in\ units\ of\ u_0$')
f_save = final_dir + '/sfigs/utot_compare_gmasks.png'
plt.savefig(f_save, dpi=512)
plt.close()
#%% NUMBER OF NEUTRAL SITES FOR EACH MASK LEVEL
m_file = final_dir + '/neutral_sites_per_mask_level.txt'
cm, ns = np.loadtxt(m_file).T
plt.figure(figsize=(7.5, 2.5))
plt.subplots_adjust(bottom=0.15, left=0.1, right=1, top=0.94)
plt.plot(cm, ns/1e8, marker='o', color='gray')
plt.xlabel('cM of neutral sites removed from telomeres')
plt.ylabel(r'number of neutral sites remaining $(\times 10^8)$')
plt.xticks(np.arange(0.1, 2.1, 0.1))
f_save = final_dir + '/sfigs/neutral_sites_per_mask.png'
plt.savefig(f_save, dpi=512)
plt.close()
#%%