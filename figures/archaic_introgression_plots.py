__author__ = 'davidmurphy'

import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir
from figures.other_code.summary_slide import cst_from_fldr

# u0 for scaling parameter estimates
u0 = 1.4e-08

#%% COMPARE BMAP FILES AND ARCHAIC PROBS DIRECTLY


def get_archaic_bvals(ch):
    cst = ChromStruct(ch)

    # create an array of b values for each position in chrom
    fb = root_dir + '/result/final_files/cadd93/bmap/{}.bmap.txt'.format(ch)
    b_arr = np.zeros(shape=cst.chlen)
    start = 0
    for (b, l) in np.loadtxt(fb):
        end = int(start + l)
        b_arr[start:end] = b
        start = end

    # get frequency archaic and mean b value for each segment in archaic map
    # fa = root_dir + '/data/snps/ArchIE-calls-YRI-MSL/YRI-freq-map.bed'
    fa = root_dir + '/data/snps/YRI.NA18933.mean.probs.txt'
    res = []
    with open(fa, 'r') as f:
        for line in f:
            line = line.split()
            c, i, j, p = line[0], int(line[1]), int(line[2]), float(line[3])
            if c == ch:
                mean_b = b_arr[i:j].mean()
                res.append((mean_b, p))

    return np.array(res)


#%% LOAD ARCHAIC INTROGRESSION LEVELS SORTED BY B


def get_archaic(fldr, pop, num):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_arch = fdir + 'predsort_archaic_{}_n{}.txt'.format(pop, num)
    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    pred, ai = np.loadtxt(f_arch).T
    pred /= pi0

    return pred, ai


def archaic_introgression_plot(fldr, pop, ttl, num):
    pred, ai = get_archaic(fldr, pop, num)
    fig_dir = root_dir + '/result/final_files/sfigs/'
    f_save = fig_dir + 'fig_S44.{}_archaic_{}_n{}.png'.format(fldr, ttl, num)
    plt.figure(figsize=(3,3))
    plt.subplots_adjust(top=0.93, right=0.99, hspace=0.1, bottom=0.14,
                        left=0.25)
    plt.title(ttl, y=0.98)
    plt.plot(pred, ai, marker='o', ms=5, lw=0, alpha=0.75, color='forestgreen')
    # plt.plot(pred, ai, alpha=0.75, color='forestgreen')

    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel('probability of archaic ancestry', labelpad=3)
    plt.xlim(0.55, 1.02)
    # plt.legend(loc='upper left', ncol=2)
    plt.savefig(f_save, dpi=512)
    plt.close()


fldr = 'cadd94_gmask_mnb_378'
# samples = ['YRI', 'MSL', 'Russian', 'Pakistan', 'YRI_indv']
# titles = ['YRI combined', 'MSL combined', 'Russian individual',
#           'Pakistani individual', 'YRI individual']
# for pop, ttl in zip(samples, titles):
#     archaic_introgression_plot(fldr, pop, ttl)
# for n in [100, 250, 500, 1000, 2000]:
#     for (pop, lab) in zip(['CEU', 'CHBS'], ['European','East Asian']):
n = 100
for pop in ['YRI', 'MSL', 'CEU', 'CHBS']:
    archaic_introgression_plot(fldr, pop, pop, n)


#%% BASIC SORT ABOVE ARCHAIC WITH MULTIPLE BIN NUMBERS


def diversity_and_archaic_combined(fldr, num, pop_id):
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    sort_file = fdir + 'basic_sort_n{}.txt'.format(num)
    div, pi, pred = np.loadtxt(sort_file).T

    rst = cst_from_fldr(fldr)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    # normalize by pi0
    pi /= pi0
    pred /= pi0
    plt.figure(figsize=(4, 7))
    plt.subplots_adjust(top=1, right=1, hspace=0.1, bottom=0.07,
                        left=0.17)

    # plot standard predicted/observed on top
    plt.subplot(211)
    plt.plot(pred, pi, marker='o', ms=5, markerfacecolor='None',
             markeredgecolor='darkorange', markeredgewidth=0.9, lw=0,
             alpha=0.75)
    xtick = np.arange(0.5, 1.01, 0.1)
    ytick = np.arange(0.5, 1.3, 0.1)
    plt.xticks(xtick, color='white')
    plt.yticks(ytick)
    # plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel(r'observed $\pi/\pi_0$', labelpad=3)

    plt.plot([0, 1], [0, 1], label=r'$y=x$', color='darkslategray',
             ls='--', alpha=0.65)
    plt.text(0.52, 1.17, r'$n_{bins}=$' + str(num))
    plt.ylim(0.5, 1.25)
    plt.xlim(0.5, 1.02)
    plt.legend(loc='lower right', ncol=3)

    # plot archaic introgression on bottom
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    f_arch = fdir + 'predsort_archaic_{}_n{}.txt'.format(pop_id, num)
    ai = np.loadtxt(f_arch)[:,1]
    plt.subplot(212)
    plt.plot(pred, ai, marker='o', ms=5, lw=0, alpha=0.75, color='forestgreen')
    # plt.plot(pred, ai, alpha=0.75, color='forestgreen')
    xtick = np.arange(0.5, 1.01, 0.1)
    plt.xticks(xtick)
    plt.xlabel(r'predicted $\pi/\pi_0$', labelpad=3)
    plt.ylabel('probability of archaic ancestry', labelpad=3)
    plt.xlim(0.5, 1.02)
    plt.xticks(xtick)
    plt.text(0.52, 0.95*max(ai), pop_id)
    sdir = root_dir + '/result/final_files/sfigs/'
    f_save = sdir + '{}_{}_n{}_archaic_sortplot.png'.format(fldr, pop_id, num)
    plt.savefig(f_save, dpi=512)
    plt.close()


for n in [100, 250, 500, 1000, 2000]:
    diversity_and_archaic_combined('cadd93_gmask', n, 'CEU')
#%%
