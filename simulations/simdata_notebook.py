__author__ = 'davidmurphy'

import os, re
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import simulations.cluster_sim_config as scfg
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import adjust_arrays, load_saved
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir


#%%
dfe = default_fitness_effects
ch = 'chr22'
dst = 'rnd_0'
smp = 2
rnd = True
fudel = 0.2
udel = scfg.mean_u * fudel
alpha = 0.05
stype = 'bs1'
itr = 0

# build dict of arguments needed for cst instance
bs1cs1_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(scfg.c_anno,),
                   bdfe=(dfe,), cdfe=(dfe,), bdir=scfg.bdir,
                   cdir=scfg.cdir, tkn=scfg.tkn, chrom=ch)
cs1_args = dict(bs_annos=(), cs_annos=(scfg.c_anno,),
                bdfe=(), cdfe=(dfe,), cdir=scfg.cdir,
                tkn=scfg.tkn, chrom=ch)
bs1_args = dict(bs_annos=(scfg.b_anno,), cs_annos=(),
                bdfe=(dfe,), cdfe=(), bdir=scfg.bdir,
                tkn=scfg.tkn, chrom=ch)

if stype == 'bs1':
    use_args = bs1_args
elif stype == 'cs1':
    use_args = cs1_args
else:
    assert stype == 'bs1cs1'
    use_args = bs1cs1_args

# initialize chromstruct and set simulation stats
cst = ChromStruct(**use_args)
cst.stat.meandiv = scfg.mean_div
cst.stat.indv = smp

# build dict of arguments needed pgm instance
pgm_args = dict(nba=cst.bnum, nca=cst.cnum, nprm=len(dfe),
                itau=scfg.init_tau, ftau=scfg.true_tau, udel=udel,
                alpha=alpha)

# initialize paramgenerator to get random input params for the simulation
pgm = scfg.ParamGenerator(**pgm_args)

# get random param row index from dst label
if dst.startswith('rnd'):
    ridx = int(dst.split('_')[-1])
    tprm = pgm.random(tau=pgm_args['ftau'])[ridx]
elif dst.startswith('ufm'):
    tprm = pgm.uniform(pattern=dst, tau=pgm_args['ftau'])
else:
    raise NameError('{} distribution does not exist'.format(dst))

params = tprm

#%%
# load all arrays
bs, cs, nu, nt, dv, pl = load_saved(cst)[1:]

# # uniform mutation rate
# nu = np.full(len(dv), cst.stat.meandiv)

#%%
# load fixed arrays used to generate predicted map from params
bs, cs, nu, _, dv, _, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)

# generate single bmap from weighted average across map matrix
uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
bwt = uvec / cst.fixed.u_fix
bsx = np.exp(np.dot(bs, bwt))

#%%
# get a B values for a selection of percentiles for each B map
pcts = [0.1, 1, 2.5, 5, 10, 25, 50, 75, 90, 99, 99.9]
bpcs = np.percentile(np.exp(bs*bwt), pcts, axis=0).T
xax = np.arange(len(pcts))

#%%
npcs = 6
wdth = 0.8 / npcs
offset = -0.5 * npcs * wdth
plt.figure(figsize=(10, 6))
ustr = r'$(\mu_{del}=7.4\cdot10^{-8})$'
plt.title('B values at select percentiles across maps ' + ustr, size=16)
for pc in bpcs:
    plt.bar(xax+offset, pc, wdth, align='edge')
    offset += wdth
plt.xticks(xax, [r'$%.1f$%%' % x for x in pcts], size=12)
plt.xlabel('percentile', size=14)
plt.ylabel('B', size=14, labelpad=20)
plt.yticks(size=14)
# plt.yscale('log')
lbl = [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)]
plt.legend(lbl)
fsave = '/Users/davidmurphy/Desktop/bmap_percentiles.png'
plt.savefig(fsave, dpi=256)
plt.close()
# plt.show()

#%%
# get a map of predictions for the given params and precalc maps
pred_map = predicted_pi(params, cst, nu, bs, cs)

#%%
rhis = plt.hist(bsx, bins=100, histtype='step', cumulative=1, normed=1)
plt.show()

#%%
# fdir = root_dir + '/result/final_files/ds_00_bound_bsx0.4_n2_sampled_mar2019/'
fdir = root_dir + '/result/final_files/ds00rnd00bth0.65/'
dr = defaultdict(list)
for f in os.listdir(fdir):
    if f.endswith('.txt'):
        fname = fdir + f
        dist = re.search('rnd_\d', fname).group()
        dr[dist].append(ChromStruct(chrom='chr1', init=fname))

for dist in dr:
    # get the list of all runs for that distribution
    rlst = dr[dist]
    # get an indexed list of all likelihoods for run list
    ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
    # get indices of the top 3 best LH runs (after sorting on LH)
    top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

    # create the plot
    plt.figure(figsize=(10, 5))
    plt.subplots_adjust(left=0.08, wspace=0.65, right=0.9)
    plt.subplot(1, 9, (1, 6))
    if 'sample' in fdir:
        plt.title('background selection sampled n=2')
    else:
        plt.title('background selection deterministic')

    x = np.arange(6)
    n = 2
    w = 0.8 / n
    s = -w / 2.0

    # get the true params for the run and plot as the first variable
    tprm = 10 ** rlst[0].stat.true_params[:-1]
    tpmf = tprm * 1e8
    plt.bar(x + s, tpmf, w, color='k', label='true')
    s += w

    # plot results from the top3 indices of runs (averaged)
    pmfs = []
    for i in top3:
        r = rlst[i]
        pmfs.append(r.uvec[0])
    pmfs = np.array(pmfs)
    mpmf = np.average(pmfs, axis=0) * 1e8

    plt.bar(x + s, mpmf, w, label='best')
    plt.ylabel(r'$\mu_{del} \cdot 10^8$')
    plt.xticks(x, [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)])
    plt.xlabel('deleterious fitness effect')
    plt.legend()

    # plot udel
    ud = sum(mpmf)  # udel
    udt = sum(tpmf)  # true udel
    sp3 = plt.subplot(197)
    sp3.yaxis.tick_right()
    sp3.yaxis.set_label_position('right')
    plt.title(r'$\mu_{del} \cdot 10^8$')
    plt.bar(0, udt, 1, color='k')
    plt.bar(1, ud, 1)
    plt.xticks([])

    # plot pi0
    plt.subplot(198)
    plt.title(r'$\pi_0$')
    pi0 = rlst[top3[0]].params[-1]
    tpi0 = rlst[top3[0]].stat.true_params[-1]
    plt.bar(0, tpi0, 1, color='k')
    plt.bar(1, pi0, 1)
    plt.xticks([])
    plt.yticks([])

    # plot clh_best - clh_true
    sp4 = plt.subplot(199)
    sp4.yaxis.tick_right()
    sp4.yaxis.set_label_position('right')
    sp4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel(r'$CLH_{avg} - CLH_{true}$')
    clh, tclh = rlst[top3[0]].stat.best_lh, rlst[top3[0]].stat.true_clh
    diff = np.exp(-1e-05 * clh) - np.exp(-1e-05 * tclh)
    col = 'green' if diff > 0 else 'red'
    plt.bar(-0.4, diff, 0.4, color=col)
    plt.ylim(0.99 * diff, 1.01 * diff)
    plt.xticks([])

    fsave = fdir + dist + '.png'
    # plt.savefig(fsave, dpi=256)
    # plt.close()
    plt.show()

