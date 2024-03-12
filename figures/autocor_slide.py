__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, default_fitness_effects, root_dir, \
    human_autosomes, RunStruct, cst_from_fldr
from scipy.stats import pearsonr, spearmanr
from figures.common_functions import format_panels


#%% FIGURE S15: autocorrelations between percent cons and cons depth
use_spc = 'ape supr vert'.split()
use_pct = '9 6 2'.split()
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/cons_autocor/'

f_fmt = pdir + '{}.percent.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 5, 8)
pct_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    pct_dict[m] = np.concatenate(vals)
#
f_pct = pdir + 'percent_correlations.txt'
with open(f_pct, 'w') as f:
    for pair in [(0, 1), (0, 2), (1, 2)]:
        vals = []
        i, j = pair
        for m in morgans:
            pear = pearsonr(pct_dict[m][:, i], pct_dict[m][:, j])[0]
            line = '{}:{} {} {}\n'.format(use_pct[i], use_pct[j], m, pear)
            f.write(line)

#
f_fmt = pdir + '{}.depth.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 4, 7)
dp_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    dp_dict[m] = np.concatenate(vals)

f_depth = pdir + 'depth_correlations.txt'
with open(f_depth, 'w') as f:
    for pair in [(0,1), (0,2), (1,2)]:
        vals = []
        i, j = pair
        for m in morgans:
            pear = pearsonr(dp_dict[m][:,i], dp_dict[m][:,j])[0]
            line = '{}:{} {} {}\n'.format(use_spc[i], use_spc[j], m, pear)
            f.write(line)


#%% FIGURE S14: PLOT of cons/depth
fig = plt.figure(figsize=(6.5, 2.16))
plt.subplots_adjust(left=0.075, wspace=0.4, top=0.97, right=0.995,
                    bottom=0.155)

colors = ['darkorange', 'darkturquoise', 'purple']
xi = np.arange(len(morgans))
ylims = (0.76, 1.02)
# depth plot
ax1 = plt.subplot(131)
format_panels(ax1)
# plt.title('A', y=0.97)
for (ii, pair) in enumerate([(0,1), (0,2), (1,2)]):
    vals = []
    i, j = pair
    for m in morgans:
        pear = pearsonr(dp_dict[m][:,i], dp_dict[m][:,j])[0]
        vals.append(pear)
    lbl = r'$\rho$' + '({},{})'.format(use_spc[i], use_spc[j])
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=7.5, lw=0, marker='o',
             color=colors[ii])

plt.ylabel(r'Pearson correlation ($\rho$)', labelpad=1)
plt.yticks(x=0.04)
plt.ylim(*ylims)
# plt.xlabel('window size ' + r'($\mathrm{log_{10}(M)}$)', labelpad=1)
plt.xticks(xi, np.log10(morgans), y=0.04)
plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
           borderaxespad=0.3, ncol=1, frameon=1,
           framealpha=0.75, facecolor='white', columnspacing=0.1,
           labelspacing=0.25)

# percent depth panel
ax2 = plt.subplot(132)
format_panels(ax2)
# plt.title('B', y=0.97)
for (ii, pair) in enumerate([(0,1), (0,2), (1,2)]):
    vals = []
    i, j = pair
    for m in morgans:
        pear = pearsonr(pct_dict[m][:,i], pct_dict[m][:,j])[0]
        # spear = spearmanr(pct_dict[m][:,i], pct_dict[m][:,j])[0]
        vals.append(pear)
    lbl = r'$\rho$' + '({}%,{}%)'.format(use_pct[i], use_pct[j])
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=7.5, lw=0, marker='o',
             color=colors[ii])

# plt.ylabel('Pearson correlation', labelpad=3)
# plt.yticks(np.arange(0.95, 1, 0.01), x=0.03)
plt.ylim(*ylims)
plt.yticks(x=0.04)
# plt.xlabel('window size ' + r'($\mathrm{log_{10}(M)}$)', labelpad=1)
# plt.xlabel('window (log10(M))', labelpad=1)
plt.xticks(xi, np.log10(morgans), y=0.04)
# plt.legend(loc='lower right', frameon=1, framealpha=1, facecolor='white')
plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
           borderaxespad=0.3, ncol=1, frameon=1,
           framealpha=0.75, facecolor='white', columnspacing=0.1,
           labelspacing=0.25)

ax3 = plt.subplot(133)
format_panels(ax3)
# plt.title('C', y=0.97)
# u gamete panel
# ugam = [2.1584456076590195, 2.020426900857384, 1.877445423267905,
#         1.7066387976109214, 1.6237627739171065, 1.6413849091560109,
#         1.6114119971511733, 1.5849677684173142]
nsites = {'fish_cons91_new': 259525076, 'fish_cons92_new':230703150,
          'fish_cons93_new': 201858377, 'fish_cons94_new': 172957936,
          'fish_cons95_new': 144127578, 'fish_cons96_new': 115350333,
          'fish_cons97_new': 86982444, 'fish_cons98_new': 61784237}
folders = ['fish_cons{}_new'.format(p) for p in range(91, 99)]
ugam = [nsites[f]*cst_from_fldr(f).stat.utot[0] for f in folders]
pc_cons = range(2, 10)
plt.bar(pc_cons, ugam[::-1], color='firebrick')
# plt.ylim(0, 2.2)
plt.xticks(pc_cons, y=0.04)
# plt.ylabel(r'$U_d$ per gamete per gen.', labelpad=1)
plt.ylabel('deleterious mutation rate\nper gamete per generation', labelpad=1)
plt.yticks(x=0.05)
plt.xlabel('% conserved', labelpad=1)

plt.text(0.08, 0.9, 'a', transform=plt.gcf().transFigure, fontweight='bold')
plt.text(0.42, 0.9, 'b', transform=plt.gcf().transFigure, fontweight='bold')
plt.text(0.758, 0.9, 'c', transform=plt.gcf().transFigure, fontweight='bold')

# udel y label
fig.add_subplot(1, 3, (1,2), frame_on=False)
plt.yticks([])
plt.xticks([])
plt.xlabel('window size ' + r'$\mathrm{(log_{10}(Morgans))}$', labelpad=12)

f_save = root_dir + '/result/final_files/sfigs/fig_S15.depth_pct_corr_and_gamete_u.png'
plt.savefig(f_save, dpi=256)
plt.close()


#%% FIGURE S17: autocorrelations between ape cons and CADD at 6%
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/cons_autocor/'
f_fmt = pdir + '{}.cadd.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 2)
pct_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    pct_dict[m] = np.concatenate(vals)


#%% FIGURE S17: PLOT autocorrelations between ape cons and CADD at 6%
plt.figure(figsize=(2.5, 2.16))
plt.subplots_adjust(left=0.2, wspace=0.27, top=0.98, right=0.88,
                    bottom=0.155)

xi = np.arange(len(morgans))
ylims = (0.76, 1.02)
# CADD/fish plot
ax1 = plt.subplot(111)
format_panels(ax1)
vals = []
for m in morgans:
    pear = pearsonr(pct_dict[m][:,0], pct_dict[m][:,1])[0]
    vals.append(pear)
lbl = r'$\rho$' + '({},{})'.format('phastCons', 'CADD')
plt.plot(xi, vals, label=lbl, alpha=0.75, ms=7.5, lw=0, marker='o', color='k')
plt.ylabel(r'Pearson correlation ($\rho$)', labelpad=1)
plt.yticks(x=0.04)
plt.ylim(*ylims)
plt.xlabel('window size ' + r'$\mathrm{(log_{10}(Morgans))}$', labelpad=1)
plt.xticks(xi, np.log10(morgans), y=0.04)
# plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
#            borderaxespad=0.3, ncol=1, frameon=1,
#            framealpha=0.75, facecolor='white', columnspacing=0.1,
#            labelspacing=0.25)
f_save = root_dir + '/result/final_files/sfigs/fish.CADD.autocorr.png'
plt.savefig(f_save, dpi=512)
plt.close()


#%% FIGURE S14: autocorrelations between ape 6% exonic and nonexonic
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/cons_autocor/'
f_fmt = pdir + '{}.exnex.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 2, 3)
an_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    an_dict[m] = np.concatenate(vals)


#%% FIGURE S14: PLOT of autocorrlations between vert 6% ex/nex
plt.figure(figsize=(2.5, 2.16))
plt.subplots_adjust(left=0.2, wspace=0.27, top=0.98, right=0.88,
                    bottom=0.155)

colors = ['darkorange', 'darkturquoise', 'purple']
labels = [r'$con_a$', r'$con_e$', r'$con_n$']

xi = np.arange(len(morgans))
ylims = (0.15, 1.02)
# ex/nex plot
ax1 = plt.subplot(111)
format_panels(ax1)
for (ii, pair) in enumerate([(0,1), (0,2), (1,2)]):
    vals = []
    i, j = pair
    for m in morgans:
        pear = pearsonr(an_dict[m][:,i], an_dict[m][:,j])[0]
        vals.append(pear)
    lbl = r'$\rho$' + '({},{})'.format(labels[i], labels[j])
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=7.5, lw=0, marker='o',
             color=colors[ii])
plt.ylabel(r'Pearson correlation ($\rho$)', labelpad=1)
plt.yticks(x=0.04)
plt.ylim(*ylims)
# plt.xlabel('window (log10(M))', labelpad=1)
plt.xlabel('window size ' + r'$\mathrm{(log_{10}(Morgans))}$', labelpad=1)

plt.xticks(xi, np.log10(morgans), y=0.04)
plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
           borderaxespad=0.3, ncol=1, frameon=1,
           framealpha=0.75, facecolor='white', columnspacing=0.1,
           labelspacing=0.25)
f_save = root_dir + '/result/final_files/sfigs/fig_S21.vert.exnex.autocorr.png'
plt.savefig(f_save, dpi=256)
plt.close()


#%% (NOT USED)
plt.figure(figsize=(2.5, 2.5))
plt.subplots_adjust(left=0.18, wspace=0.27, top=1, right=1, bottom=0.15)

xi = np.arange(len(morgans))
labels = [r'$cons_a$', r'$cons_e$', r'$cons_n$']
for pair in [(0,1), (0,2), (1,2)]:
    vals = []
    i, j = pair
    for m in morgans:
        pear = pearsonr(an_dict[m][:,i], an_dict[m][:,j])[0]
        vals.append(pear)
    lbl = r'$\rho$' + '({},{})'.format(labels[i], labels[j])
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=7.5, lw=0, marker='o')

plt.ylabel('Pearson correlation', labelpad=3)
plt.yticks(x=0.03)
# plt.yticks(np.arange(0.90, 1.01, 0.02), x=0.03)
plt.ylim(0.40, 1.05)
plt.xlabel('Window size (log10(Morgans))', labelpad=3)
plt.xticks(xi, np.log10(morgans), y=0.03)
plt.legend(loc='lower right', frameon=1, framealpha=1, facecolor='white')
f_save = root_dir + '/result/final_files/sfigs/fig_S21.vert.exnex.autocorr.png'
plt.savefig(f_save, dpi=256)
plt.close()


#%% autocorrelations between 7% CADD and 7% CADD marked HC-runsubs
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/hcsubs_cons_vs_cadd_autocor/'
f_fmt = pdir + '{}.cscs.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 2)
sb_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    sb_dict[m] = np.concatenate(vals)

pdir = root_dir + '/data/phast/hcsubs_cons_autocor/'
f_fmt = pdir + '{}.bscs.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 2)
cn_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    cn_dict[m] = np.concatenate(vals)

pdir = root_dir + '/data/phast/hcsubs_cadd_autocor/'
f_fmt = pdir + '{}.bscs.gmap.dist.bin_{:.2e}.txt'
use_cols = (1, 2)
cd_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f, usecols=use_cols)
        vals.append(arr)

    cd_dict[m] = np.concatenate(vals)

#%% PLOT of autocorrelations between 7% CADD and 7% CADD marked HC-subs
plt.figure(figsize=(3.5, 3.5))
plt.subplots_adjust(left=0.2, wspace=0.27, top=1, right=1, bottom=0.15)
xi = np.arange(len(morgans))

# Ape cons anno to ape cons subs correlation
vals = []
labs = ['ape', 'ape-subs']
for m in morgans:
    pear = pearsonr(cn_dict[m][:,0], cn_dict[m][:,1])[0]
    vals.append(pear)
lbl = r'$\rho$' + '({},{})'.format(*labs)
plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o')

# CADD anno to CADD subs correlation
vals = []
labs = ['CADD', 'CADD-subs']
for m in morgans:
    pear = pearsonr(cd_dict[m][:,0], cd_dict[m][:,1])[0]
    vals.append(pear)
lbl = r'$\rho$' + '({},{})'.format(*labs)
plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o')

# subs correlation
vals = []
labs = ['CADD-subs', 'ape-subs']
for m in morgans:
    pear = pearsonr(sb_dict[m][:,0], sb_dict[m][:,1])[0]
    vals.append(pear)
lbl = r'$\rho$' + '({},{})'.format(*labs)
plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o')

plt.ylabel('Pearson correlation')
plt.xlabel('log10(M)')
plt.xticks(xi, np.log10(morgans))
plt.ylim(0.645, 0.991)
plt.legend(loc='lower right', frameon=1, framealpha=1, facecolor='white')
# f_save = root_dir + '/data/phast/CADD7pct.bscs.autocorr.png'
# f_save = root_dir + '/data/phast/ape7pct.bscs.autocorr.png'
f_save = root_dir + '/data/phast/BS_vs_CS_annos.autocorr.png'

plt.savefig(f_save, dpi=256)
plt.close()

#%% FIGURE S19: correlations between CADD 7% segments and CADD 2-10% HC subs
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/hcsubs_cons_autocor/'
f_fmt = pdir + '{}.bscs.gmap.dist.bin_{:.2e}.txt'
sb_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f)[:,1:]
        vals.append(arr)

    sb_dict[m] = np.concatenate(vals)

#%% FIGURE S19: PLOT correlations between CADD 7% segments and CADD 2-10% HC subs
plt.figure(figsize=(2.5, 2.16))
plt.subplots_adjust(left=0.2, wspace=0.27, top=0.98, right=0.975,
                    bottom=0.155)
ax1 = plt.subplot(111)
format_panels(ax1)
labels = [r'$cons$', r'$cons_{ex}$', r'$cons_{nex}$']
ylims = (0.65, 1.02)

xi = np.arange(len(morgans))
pct = range(91, 99)
# pct = [90, 98]
labs = ['{:>3}%'.format(100-p) for p in pct]
# cols = ['orange','purple']
cols = ['mediumpurple','mediumvioletred',  'red', 'gold', 'rosybrown',
        'maroon', 'turquoise', 'firebrick', 'dodgerblue']
for i in [1, 8]:
    vals = []
    for m in morgans:
        pear = pearsonr(sb_dict[m][:,0], sb_dict[m][:,i])[0]
        vals.append(pear)
    lbl = labs[i-1]
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o',
             color=cols[i-1])

plt.ylabel(r'Pearson correlation ($\rho$)', labelpad=1)
plt.yticks(x=0.04)
plt.ylim(*ylims)
plt.xlabel('window size ' + r'$\mathrm{(log_{10}(Morgans))}$', labelpad=1)
plt.xticks(xi, np.log10(morgans), y=0.04)
plt.legend(loc='lower right', handletextpad=-0.3, borderpad=0.2,
           borderaxespad=0.3, ncol=1, frameon=1,
           framealpha=0.75, facecolor='white', columnspacing=0.1,
           labelspacing=0.25)
f_save = root_dir + '/result/final_files/sfigs/fish_cons_hc_subs.autocorr.png'
plt.savefig(f_save, dpi=512)
plt.close()


#%% correlations between CADD 7% segments and ape 2-10% HC subs
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/caddBS_apeCS_all_pct_autocor/'
f_fmt = pdir + '{}.bscs.gmap.dist.bin_{:.2e}.txt'
sb_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f)[:,1:]
        vals.append(arr)

    sb_dict[m] = np.concatenate(vals)

#%% PLOT correlations between CADD 7% segments and ape 2-10% HC subs
plt.figure(figsize=(3.5, 3.5))
plt.subplots_adjust(left=0.2, wspace=0.27, top=0.9, right=0.95, bottom=0.15)
plt.title('CADD 7% segments vs. ape 2-10% subs', fontsize=11)
xi = np.arange(len(morgans))
pct = range(90, 99)
labs = ['{}%'.format(100-p) for p in pct]
cols = ['mediumvioletred', 'red', 'orange', 'gold', 'rosybrown',
        'maroon', 'turquoise', 'dodgerblue', 'purple']
# depth plot
for i in range(1, 10):
    vals = []
    for m in morgans:
        pear = pearsonr(sb_dict[m][:,0], sb_dict[m][:,i])[0]
        vals.append(pear)
    lbl = labs[i-1]
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o',
             color=cols[i-1])

plt.ylabel('Pearson correlation')
# plt.yticks(np.arange(0.90, 1.01, 0.02))
plt.ylim(0.55,0.99)
plt.xlabel('log10(M)')
plt.xticks(xi, np.log10(morgans))
plt.legend(loc='lower right', frameon=1, framealpha=0.75, facecolor='white',
           ncol=3)
f_save = root_dir + '/data/phast/caddBS.apeCS.autocorr.png'
plt.savefig(f_save, dpi=256)
plt.close()

#%% correlations between ape 7% segments and ape 2-10% HC subs
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/apeBS_apeCS_all_pct_autocor/'
f_fmt = pdir + '{}.bscs.gmap.dist.bin_{:.2e}.txt'
sb_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f)[:,1:]
        vals.append(arr)

    sb_dict[m] = np.concatenate(vals)

#%% PLOT correlations between ape 7% segments and ape 2-10% HC subs
plt.figure(figsize=(3.5, 3.5))
plt.subplots_adjust(left=0.2, wspace=0.27, top=0.9, right=0.95, bottom=0.15)
plt.title('Ape 7% segments vs. ape 2-10% subs', fontsize=11)
xi = np.arange(len(morgans))
pct = range(90, 99)
labs = ['{}%'.format(100-p) for p in pct]
cols = ['mediumvioletred', 'red', 'orange', 'gold', 'rosybrown',
        'maroon', 'turquoise', 'dodgerblue', 'purple']
# depth plot
for i in range(1, 10):
    vals = []
    for m in morgans:
        pear = pearsonr(sb_dict[m][:,0], sb_dict[m][:,i])[0]
        vals.append(pear)
    lbl = labs[i-1]
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o',
             color=cols[i-1])

plt.ylabel('Pearson correlation')
# plt.yticks(np.arange(0.90, 1.01, 0.02))
plt.ylim(0.55,0.99)
plt.xlabel('log10(M)')
plt.xticks(xi, np.log10(morgans))
plt.legend(loc='lower right', frameon=1, framealpha=0.75, facecolor='white',
           ncol=3)
f_save = root_dir + '/data/phast/apeBS.apeCS.autocorr.png'
plt.savefig(f_save, dpi=256)
plt.close()

#%% correlations between ape 7% segments and CADD 2-10% HC subs
morgans = [1e-4, 1e-3, 1e-2, 1e-1]
pdir = root_dir + '/data/phast/apeBS_caddCS_all_pct_autocor/'
f_fmt = pdir + '{}.bscs.gmap.dist.bin_{:.2e}.txt'
sb_dict = {}
for m in morgans:
    vals = []
    for ch in human_autosomes:
        f = f_fmt.format(ch, m)
        arr = np.loadtxt(f)[:,1:]
        vals.append(arr)

    sb_dict[m] = np.concatenate(vals)

#%% PLOT correlations between ape 7% segments and CADD 2-10% HC subs
plt.figure(figsize=(3.5, 3.5))
plt.subplots_adjust(left=0.2, wspace=0.27, top=0.9, right=0.95, bottom=0.15)
plt.title('Ape 7% segments vs. CADD 2-10% subs', fontsize=11)
xi = np.arange(len(morgans))
pct = range(90, 99)
labs = ['{}%'.format(100-p) for p in pct]
cols = ['mediumvioletred', 'red', 'orange', 'gold', 'rosybrown',
        'maroon', 'turquoise', 'dodgerblue', 'purple']
# depth plot
for i in range(1, 10):
    vals = []
    for m in morgans:
        pear = pearsonr(sb_dict[m][:,0], sb_dict[m][:,i])[0]
        vals.append(pear)
    lbl = labs[i-1]
    plt.plot(xi, vals, label=lbl, alpha=0.75, ms=10, lw=0, marker='o',
             color=cols[i-1])

plt.ylabel('Pearson correlation')
# plt.yticks(np.arange(0.90, 1.01, 0.02))
plt.ylim(0.55, 0.99)
plt.xlabel('log10(M)')
plt.xticks(xi, np.log10(morgans))
plt.legend(loc='lower right', frameon=1, framealpha=0.75, facecolor='white',
           ncol=3)
f_save = root_dir + '/data/phast/apeBS.caddCS.autocorr.png'
plt.savefig(f_save, dpi=256)
plt.close()
#%% plot HC sub counts for ape/CADD 2-10%
f_counts = '/Users/davidmurphy/Desktop/annotated_sub_counts.txt'
perc, apes, cads = [], [], []
with open(f_counts, 'r') as f:
    for line in f:
        line = line.split()
        pct = 100 - int(line[0][:2])
        ape = float(line[1].split('=')[1])
        cad = float(line[2].split('=')[1])
        perc.append(pct)
        apes.append(ape/1e6)
        cads.append(cad/1e6)

w = 0.4
apes.reverse()
cads.reverse()
perc.reverse()
xi = np.arange(len(perc))
plt.figure(figsize=(7, 3.5))
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.99, top=0.99)
plt.bar(xi-0.4, cads, width=0.4, label='CADD', color='darkorange')
plt.bar(xi, apes, width=0.4, label='ape', color='dodgerblue')
plt.xticks(xi, ['{}%'.format(p) for p in perc])
plt.xlabel('score threshold')
plt.ylabel('HC substitution count (millions)')
plt.legend()
f_save = root_dir + '/data/phast/HCsub.counts.ape.CADD.png'
plt.savefig(f_save, dpi=256)
plt.close()

#%% calculate overlap between ape/CADD annotated HC subs for 2-10%
cdir = root_dir + '/data/csanno'
ffmt = cdir + '/{a}/{c}.{a}.npz'
avals, cvals, bothvals = [], [], []
pc = 98
for pc in range(90, 99):
    a_an = 'ape_cons{}_clean_hc_subs'.format(pc)
    c_an = 'cadd{}_hc_subs'.format(pc)
    a_only, c_only, ac = 0, 0, 0
    for c in range(1, 23):
        ch = 'chr{}'.format(c)
        a_f = ffmt.format(a=a_an, c=ch)
        c_f = ffmt.format(a=c_an, c=ch)
        ape = np.load(a_f)['pos']
        cad = np.load(c_f)['pos']
        cadinape = np.in1d(ape, cad)
        ac += np.sum(cadinape) / 1e6
        a_only += np.sum(~cadinape) / 1e6
        c_only += (len(cad) - np.sum(cadinape)) / 1e6

    avals.append(a_only)
    cvals.append(c_only)
    bothvals.append(ac)

avals.reverse()
cvals.reverse()
bothvals.reverse()

#%% PLOT overlap between CADD, ape and total counts
aheight = [a-b for a,b in zip(avals, bothvals)]
cheight = [c-b for c,b in zip(cvals, bothvals)]

xi = np.arange(len(bothvals))
plt.figure(figsize=(7, 3.5))
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.99, top=0.99)
plt.bar(xi-0.2, cvals, width=0.4, label='CADD only', color='darkorange')
plt.bar(xi+0.2, avals, width=0.4, label='ape only', color='dodgerblue')
plt.bar(xi, bothvals, width=0.8, label='shared CADD & ape',
        color='darkslategray', alpha=0.75, hatch='//')
plt.xticks(xi, ['{}%'.format(p) for p in perc])
plt.xlabel('score threshold')
plt.ylabel('HC substitution count (millions)')
plt.legend()
f_save = root_dir + '/data/phast/HCsub.counts.ape.CADD.png'
plt.savefig(f_save, dpi=256)
plt.close()
#%%
pct = range(90, 99)
for p in pct:
    cst = cst_from_fldr('cadd{}_BSCS'.format(p))
    cst.stat.calc_stats(cst)
    atot = cst.stat.atot
    print 'CADD {}% = {}'.format(p, atot)

    cst = cst_from_fldr('ape_cons{}_clean_BSCS'.format(p))
    cst.stat.calc_stats(cst)
    atot = cst.stat.atot
    print 'ape {}% = {}'.format(p, atot)
#%%
cst = cst_from_fldr('ape_cons90_clean_CS_only')
cst.stat.calc_stats(cst)
print cst.stat.atot
#%%
fldr = root_dir +  '/result/final_files/ape_cons90_clean_BSCS/'
finits = [fldr + f for f in os.listdir(fldr) if f.startswith('YRI.')]
rlst = [ChromStruct('chr1', init=f) for f in finits]
# assign indices to each runstruct
ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
# get indices of the top 3 best LH runs (after sorting on LH)
top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

#%%
fldr_1 = 'ape_cons90_4pt'
fldr_2 = 'ape_cons90_BSCS_4pt'
fldr_3 = 'ape_cons90_BS_CSNS_4pt'
f_save = root_dir + '/result/final_files/{}/compare_3.png'.format(fldr_1)
xtck = [-3.5, -3.0, -2.5, -2.0]
xi = np.arange(4)
xlim = (-0.5, 3.5)

cst_1 = cst_from_fldr(fldr_1)
cst_2 = cst_from_fldr(fldr_2)
cst_3 = cst_from_fldr(fldr_3)
r1 = np.loadtxt(root_dir + '/result/final_files/{}/rsq.log'.format(fldr_1))
r2 = np.loadtxt(root_dir + '/result/final_files/{}/rsq.log'.format(fldr_2))
r3 = np.loadtxt(root_dir + '/result/final_files/{}/rsq.log'.format(fldr_3))

# plot(2,8,(1,3)); (284), (285), (286), (2,4,4)

plt.figure(figsize=(10, 3))
plt.subplots_adjust(left=0.07, wspace=1.5, right=1, bottom=0.15, top=0.98,
                    hspace=0.75)

# BS params
plt.subplot(2,8,(1,3))
# plt.title('BS parameters')
plt.bar(xi-0.3, cst_1.stat.uvec[0] / 1.4e-8, 0.3, label='BS only', color='darkorange')
plt.bar(xi, cst_2.stat.uvec[0] / 1.4e-8, 0.3, label='BS+CS: NS ape 10%', color='dodgerblue')
plt.bar(xi+0.3, cst_3.stat.uvec[0] / 1.4e-8, 0.3, label='BS+CS: NS all', color='purple')

plt.ylim(0., 0.9)
plt.ylabel(r'$\mu_{del}/\mu_{tot}$')
plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck])
plt.xlim(xlim)
plt.xlabel('deleterious fitness effect')
plt.legend(ncol=1)

# CS params
csaxis = plt.subplot(2,8,(9,11))
# plt.title('CS parameters')
plt.bar(xi, cst_2.stat.avec[0], 0.3 , color='dodgerblue')
plt.bar(xi+0.3, cst_3.stat.avec[0], 0.3, color='purple')

# plt.bar(xi+0.2, cst_2.stat.avec[1], 0.4, label='CS: other')
plt.ylabel(r'$\alpha$')
plt.xlabel('adaptive fitness effect')
plt.xticks(xi, [r'$10^{%.1f}$' % x for x in xtck])
plt.xlim(xlim)
csaxis.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# plt.xticks(xi, ['' for _ in xtck])
plt.legend(loc='upper left')

plt.subplot(184)
plt.ylabel('mutation rate (in units of ' + r'$\mu_{tot}$' + ')')
plt.bar(0, cst_1.stat.uvec[0].sum() / 1.4e-8, 1, color='darkorange')
plt.bar(1, cst_2.stat.uvec[0].sum() / 1.4e-8, 1, color='dodgerblue')
plt.bar(2, cst_3.stat.uvec[0].sum() / 1.4e-8, 1, color='purple')

plt.xticks(color='none')

mpi01 = cst_1.params[-1] / cst_1.fixed.tau_init
mpi02 = cst_2.params[-1] / cst_2.fixed.tau_init
mpi03 = cst_3.params[-1] / cst_3.fixed.tau_init

plt.subplot(185)
plt.ylabel('average reduction in heterozygosity')
plt.bar(0, mpi01, 1, label='BS only', color='darkorange')
plt.bar(1, mpi02, 1, label='BS+CS', color='dodgerblue')
plt.bar(2, mpi03, 1, label='BS+CS', color='purple')

plt.xticks(color='none')

plt.subplot(186)
plt.ylabel(r'$- \Delta CLLH\ (x10^5)$')
clh_list = [cst_1.stat.best_lh, cst_2.stat.best_lh, cst_3.stat.best_lh]
clh_list = [1e5 * (cl / min(clh_list) - 1) for cl in clh_list]
plt.bar(0, clh_list[0], 1, label='BS only', color='darkorange')
plt.bar(1, clh_list[1], 1, label='BS+CS', color='dodgerblue')
plt.bar(2, clh_list[2], 1, label='BS+CS', color='purple')

plt.xticks(color='none')

plt.subplot(1,8,(7,8))
plt.ylabel('variance explained ' + r'$(R^2)$')
rsq_list = [r1, r2, r3]
ws = np.log10(rsq_list[0][:,0])
wmsk = (ws < 6.4)
ws = ws[wmsk]
colors = ['darkorange', 'dodgerblue', 'purple']
label_list = ['BS only', 'BS+CS', 'BS+CS']
for (i, rs) in enumerate(rsq_list):
    # save the 1Mb rsq vals
    rsq = rs[:,1][wmsk]
    plt.plot(ws, rsq, color=colors[i], marker='o',
             linewidth=0, alpha=0.75, label=label_list[i])
plt.xlabel('window size (bp)')
# plot subset of x ticks
xtck = [4, 5, 6]
plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4,5,6]])

plt.savefig(f_save, dpi=256)
plt.close()

#%%
cdir = root_dir + '/data/csanno'
ffmt = cdir + '/{a}/chr{c}.{a}.npz'
acon = 'nonsyn_cadd90_hc_subs'
aall = 'nonsyn'
ncon = 0
nall = 0
for c in range(1, 23):
    ncon += np.load(ffmt.format(a=acon, c=c))['pos'].size
    nall += np.load(ffmt.format(a=aall, c=c))['pos'].size
print 'all={} cons={}'.format(nall, ncon)
#%%