__author__ = 'davidmurphy'

import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from classes.runstruct import root_dir, RunStruct, human_autosomes


def chr_map(rdir):
    fdir = root_dir + '/result/final_files/{}/'.format(rdir)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    clist = [f for f in os.listdir(fdir) if ('chr1' in f) and
             (f.endswith('.txt'))]
    assert len(flist) == len(clist) == 1
    f_init = fdir + flist[0]
    f_name = fdir + clist[0]
    rst = RunStruct(init=f_init)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    f_save = f_name.replace('.txt', '.png')
    prd, obs, num = np.loadtxt(f_name).T
    xi = np.arange(0, prd.size / 2.0, 0.5)
    if 'cthresh' in rdir:
        pi0 = prd[-1]
    pi_mean = np.nanmean(obs)
    prd /= pi_mean
    obs /= pi_mean
    obs[(obs < 0.25) | (obs > 1.75)] = np.nan
    plt.figure(figsize=(10, 2.5))
    # plt.title(rdir.split('_')[0])
    plt.title(rdir)

    plt.subplots_adjust(left=0.08, bottom=0.18, right=0.98)
    plt.plot(xi, obs, label='observed', color='darkslategray', lw=1.5)
    plt.plot(xi, prd, label='predicted', color='fuchsia', lw=1.5, alpha=0.8)
    plt.ylabel(r'$\pi/\bar{\pi}$')
    plt.ylim(0.1,1.75)
    plt.xlabel('chr1 position (Mb)')
    plt.xticks(range(25, 250, 50))
    plt.xlim(0, xi[-1])
    plt.legend()
    plt.savefig(f_save, dpi=256)
    plt.close()


def sort_map(rdir):
    fdir = root_dir + '/result/final_files/{}/'.format(rdir)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    slist = [f for f in os.listdir(fdir) if ('predsort' in f) and
             (f.endswith('.txt'))]
    assert len(flist) == len(slist) == 1
    f_init = fdir + flist[0]
    f_name = fdir + slist[0]
    rst = RunStruct(init=f_init)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    f_in = f_name
    f_save = f_in.replace('.txt', '.png')
    div, pi, pred = np.loadtxt(f_in).T
    xi = np.arange(div.size)
    if 'cthresh' in rdir:
        pi0 = pred[-1]
    pi /= pi0
    pred /= pi0
    plt.figure(figsize=(4,4))
    plt.subplots_adjust(left=0.15, bottom=0.15)
    # plt.title(rdir.split('_')[0])
    plt.title(rdir)
    plt.axhline(y=1, color='k', ls='--', alpha=0.8)
    plt.text(50, 1.03, 'w/o background selection', ha='center', va='center')
    meanpi = pi.mean()
    plt.axhline(y=meanpi, color='darkslategray', ls=':', alpha=0.8)
    plt.text(60, meanpi-0.03, 'mean diversity', ha='left', va='center',
             color='darkslategray')

    plt.plot(pred, pi, label='observed', color='darkslategray')
    plt.plot(pred, pred, label='predicted', color='fuchsia')
    plt.ylabel(r'$\pi/\pi_0$', fontsize=14)
    plt.ylim(0.5, 1.2)
    plt.xlim(0.5, 1.02)
    plt.xlabel('strength of background selection')
    plt.savefig(f_save, dpi=256)
    plt.close()


def rsq_map(rdir):
    f_mcv = '{}/result/rsq_maps/old_data_rsq/phylo-2/' \
           'mcvicker_map_BS1_CS0_161129211053_rsqmap.txt'.format(root_dir)
    f_ape = root_dir + '/result/final_files/{}/rsq.log'.format(rdir)
    f_exonly = root_dir + '/result/final_files/excons_only/rsq.log'
    f_save = f_ape.replace('rsq.log', 'rsq.png')
    w, mcv = np.loadtxt(f_mcv)[:16].T
    ape = np.loadtxt(f_ape)[:16,1]
    ex = np.loadtxt(f_exonly)[:16,1]

    w = np.log10(w)
    plt.figure(figsize=(5,5))
    plt.plot(w, mcv, label='McVicker', marker='o', lw=0, color='cyan')
    plt.plot(w, ape, label='ape 5% conserved', marker='o', lw=0, color='fuchsia')
    plt.plot(w, ex, label='ape 5% conserved exonic', marker='o', lw=0, color='purple')

    plt.xlabel('window size')
    xtck = [4, 5, 6]
    plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4,5,6]])
    plt.title('variance explained')
    plt.ylabel(r'$R^2$')
    plt.legend(prop=dict(size=9), loc='upper left')
    plt.savefig(f_save, dpi=256)
    plt.close()


def collate_plot(rdir):
    fdir = root_dir + '/result/final_files/{}/'.format(rdir)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    f_init = fdir + flist[0]
    rst = RunStruct(init=f_init)
    pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]

    f_ape = root_dir + '/result/final_files/{}/nonsyn.collate.5.00e-03width.npy'.format(rdir)
    f_save = f_ape.replace('.npy', '.png')
    f_mcv = root_dir + '/result/collate/YRI.mcvicker.ref.clean.BS1.1.CS0.0.nonsyn.collated.npy'

    bins, div, pi, pr, cnts = np.load(f_ape).T
    _, _, _, mcv, mcv_cnts = np.load(f_mcv).T
    # bins, pi, div, pred, cnts = np.load(f_in).T
    obs = pi / cnts
    prd = pr / cnts
    mcv /= mcv_cnts
    # y_max = pi0
    y_max = obs.mean()
    obs /= y_max
    prd /= y_max
    mcv /= y_max
    plt.figure(figsize=(5, 3))
    plt.title(rdir)
    plt.subplots_adjust(bottom=0.15)
    plt.plot(bins, obs, label='observed', color='darkslategray')
    plt.plot(bins, prd, label='predicted', color='fuchsia')
    # plt.plot(bins, mcv, label='McVicker', color='cyan')
    plt.xlabel('distance to nearest NS substitution (cM)')
    # plt.ylabel(r'$\pi/\pi_0$', fontsize=14)
    plt.ylabel(r'$\pi/\bar{\pi}$', fontsize=14)

    plt.legend()
    plt.savefig(f_save, dpi=256)
    plt.close()


#%%
spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()
pct = range(90, 97, 1)
for pc in pct:
    rdir = 'ape_cons{}_clean'.format(pc)
    # rdir = 'ape_cons95_clean'
    chr_map(rdir)
    sort_map(rdir)
    collate_plot(rdir)

#%%
spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()
flist = ['{}_cons95_clean'.format(sp) for sp in spec]
rsq_vals = {}
rsq_rank = dict((sp, 0) for sp in spec)
for sp in spec:
    fldr = '{}_cons95_clean'.format(sp)
    fdir = root_dir + '/result/final_files/{}/'.format(fldr)
    rsq_vals[sp] = np.loadtxt(fdir+'rsq.log')[:16,1]

ranks = []
for i in range(len(rsq_vals['fish'])):
    ri = [rsq_vals[sp][i] for sp in spec]
    si = np.argsort(ri)
    ranks.append([spec[j] for j in si])
    ri_max = max(ri)
    for sp in spec:
        if rsq_vals[sp][i] == ri_max:
            rsq_rank[sp] += 1

plt.figure(figsize=(5,5))
plt.subplots_adjust(bottom=0.25)
plt.bar(range(7), [rsq_rank[sp] / 16.0 for sp in spec])
plt.ylabel('% '+r'$R^2$' + ' rank=1')
plt.xticks(range(7), spec, rotation=90)
f_save = root_dir + '/result/final_files/rank_rsq_clean.png'
plt.savefig(f_save, dpi=256)
plt.close()
# plt.show()

 #%%
# rdir = 'cthresh_001'
# llist = '20 30 40 50 60 65 70 75 80 85'.split()
# flist = ['nothresh_4prm']
# flist += ['ape95_bth{}'.format(l) for l in llist]
mbsuf = '0001 001 01 05 1 15 20 25'.split()
ctsuf = '001 01 05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100'
# flist = ['ape95_minbs_{}'.format(sf) for sf in mbsuf]
flist = ['cthresh_{}'.format(l) for l in ctsuf.split()]

for rdir in flist:
    chr_map(rdir)
    sort_map(rdir)



#%% SORT PLOT MCVICKER
f_mcv = root_dir + '/result/sort_map/alternateclean/YRI.mcvicker.ref.' \
                   'clean.BS1.1.CS0.0.sorted.100bins.npy'
f_save = root_dir + '/result/final_files/mcvicker_sorted.png'
# div, pi, pred = np.loadtxt(f_thresh).T
div, pi, pred = np.load(f_mcv).T
xi = np.arange(div.size)
pi /= pred.max()
pred /= pred.max()
plt.figure(figsize=(4,4))
plt.subplots_adjust(left=0.15, bottom=0.15)
plt.axhline(y=1, color='k', ls='--', alpha=0.8)
# plt.text((0.5 * num_bins), 1.005, 'w/o background selection', ha='center', va='bottom', fontsize=18)
plt.text(50, 1.03, 'w/o background selection', ha='center', va='center')
meanpi = 0.82
plt.axhline(y=meanpi, color='darkslategray', ls=':', alpha=0.8)
plt.text(60, 0.79, 'mean diversity', ha='left', va='center', color='darkslategray')

plt.plot(xi, pi, label='observed', color='darkslategray')
plt.plot(xi, pred, label='McVicker', color='darkcyan')
plt.ylabel('scaled diversity')
plt.xlabel('strength of background selection')
plt.savefig(f_save, dpi=256)
plt.close()


#%% DFE PLOT
fdir = root_dir + '/result/final_files/'
fldrs = [fdir + 'ape95_euarchontoglires35_filtered/']
fmt = 'ape_cons{}_euarchontoglires_neut35_filtered/'
for pc in [94, 93, 92, 91]:
    fldrs.append(fdir + fmt.format(pc))

for (i, fl) in enumerate(fldrs):
    flist = [f for f in os.listdir(fl) if f.endswith('composite.txt')]
    assert len(flist) == 1
    f_init = fl + flist[0]
    rst = RunStruct(f_init)
    uvec = rst.stat.uvec[0]

plt.figure(figsize=(5, 3))


#%% U DELETERIOUS PLOT
col = ['DarkTurquoise', 'Orange', 'DarkOrange', 'Coral', 'Tomato', 'OrangeRed']
fig = plt.figure(figsize=(10, 7))
fig.subplots_adjust(left=0.14, bottom=0.19, right=1, top=0.95, wspace=0,
                    hspace=0)
baseline = r'$\mathrm{Kong\ (2012)}$'
totalrate = r'$\mathrm{Total\ rate\ (?)}$'

plt.bar(left=0, height=7.4, width=0.9, color='DarkTurquoise', align='edge',
        label='McVicker\'s estimate: conserved exonic')
plt.bar(left=1, height=6.21, width=0.9, color='Purple', align='edge',
        label='Our estimate: conserved exonic')
plt.bar(left=2, height=1.07, width=0.9, color='Fuchsia', align='edge',
        label='Our estimate: all conserved')
plt.bar(left=3, height=1.03, width=0.9, alpha=0.8, color=col[2], align='edge')
plt.bar(left=4, height=0.876, width=0.9, alpha=0.8, color=col[3], align='edge')
plt.bar(left=5, height=0.788, width=0.9, alpha=0.8, color=col[4], align='edge')
plt.bar(left=6, height=0.71, width=0.9, alpha=0.8, color=col[5], align='edge')
plt.axhline(y=1.2, xmin=-0.1, xmax=2, ls='--', lw=2, color='DarkSlateGray')
plt.axhline(y=2.4, xmin=-0.1, xmax=2, ls=':', lw=2, color='Gray', alpha=0.75)
plt.text(0.45, 7.5, r'$\mathrm{McVicker}$', fontsize=22, ha='center',
         va='bottom')
plt.text(1.45, 6.31, r'$\mathrm{5\%}$' '\n' r'$\mathrm{excon}$', fontsize=22,
         ha='center', va='bottom')
plt.text(4, 1.25, 'Kong (2012)', fontsize=24, ha='left', va='bottom',
         color='darkslategray')
plt.text(4, 2.45, 'total rate (?)', fontsize=24, ha='left', va='bottom',
         color='DimGray')

plt.annotate(r'$\mathrm{\%\ conserved}$', xycoords='figure fraction',
             xy=(0.68, 0.1), xytext=(0.68, 0.01), fontsize=24, ha='center',
             va='bottom',
             arrowprops=dict(arrowstyle='-[, widthB=6.55, lengthB=0.5', lw=2.0))
plt.ylabel(
    r'$\mathrm{deleterious\ mutation\ rate}$' '\n' r'$\mathrm{per\ bp/gen\ (x10^{-8})}$',
    fontsize=24, labelpad=15)

# plt.yticks(np.arange(1, 9), ['{:.1f}'.format(y) for y in range(1, 9)], **tickprop)
plt.yticks(np.arange(1, 9), fontsize=22)
plt.yticks(fontsize=24)
plt.xticks([])
plt.ylim(0, 8.5)
# plt.xlabel(r'$\mathrm{fraction\ conserved}$', fontsize=32, labelpad=25)
plt.xticks(np.arange(0.5, 7, 1), ['', '', '5%', '6%', '7%', '8%', '9%'],
           fontsize=24)
f_save = root_dir + '/result/final_files/udel_vs_anno.png'
plt.savefig(f_save, dpi=256)
plt.close()

#%% EX/NEX DFE PLOT
f_init = root_dir + '/result/final_files/ape_exnex/YRI.ape_exnex.BS2.4.' \
                    'CS0.0.iprm_05.bth_0.650.191001143948.final.composite.txt'
rst = RunStruct(init=f_init)
ex, nex = rst.uvec
xi = np.arange(4)
plt.figure(figsize=(4,4))
plt.subplot(211)
plt.bar(xi, ex*1e8)
plt.subplot(212)
plt.bar(xi, nex*1e8)
plt.show()

#%% R SQAURED COMPARE CS-ONLY RESULTS
f_ape = root_dir + '/result/final_files/ape95_euarchontoglires35_filtered/rsq.log'
f_csns = root_dir + '/result/final_files/cs_only_nonsyn/rsq.log'
f_cscons = root_dir + '/result/final_files/hc_derived_cons_only/rsq.log'
f_save = root_dir + '/result/final_files/compare_CS_only_rsq.png'
w, ape = np.loadtxt(f_ape)[:16].T
ns = np.loadtxt(f_csns)[:16,1]
hc = np.loadtxt(f_cscons)[:16,1]

w = np.log10(w)
plt.figure(figsize=(5,5))
plt.plot(w, ape, label='ape 5% conserved', marker='o', lw=0, color='fuchsia')
plt.plot(w, hc, label='CS HC-derived ape 5% cons', marker='o', lw=0, color='purple')
plt.plot(w, ns, label='CS nonsyn', marker='o', lw=0, color='cyan')

plt.xlabel('window size')
xtck = [4, 5, 6]
plt.xticks(xtck, [r'$10^{%.1f}$' % x for x in [4,5,6]])
plt.title('variance explained')
plt.ylabel(r'$R^2$')
plt.legend(prop=dict(size=9), loc='upper left')
plt.savefig(f_save, dpi=256)
plt.close()

#%%
# for pc in [91, 92, 93, 94]:
# spec = 'ape fish primate prosimian euarchontoglires laurasiatheria mammal'
# for cons in spec.split():
cons = 'ape'
pc = 95
npct = 0.35
# cons = 'ape'
ncons = 'euarchontoglires'
an_1 = '{}_cons{}_exonic'.format(cons, pc)
an_2 = '{}_cons{}_nonexonic'.format(cons, pc)
# tk ='ape_exnex'
tk = an_1
sl = False
wn = 7000
ns = 8
bs_tkn = tk
cs_tkn = 'hc_derived_cons'
# for neut in 'LWK CEU JPT'.split():
neut = 'YRI'
dpop = None
dfe_6pt = np.power(10.0, np.arange(-4.5, -1.5, 0.5))
a = RunStruct(nspc=ns, wind=wn, slid=sl, tkn=tk, ncon=ncons, npct=npct,
              cons=cons, neut=neut, dpop=dpop,
              bs_annos=(an_1,),
              # cs_annos=(cs_tkn,),
              bdir='ape_exnex',
              # cdir=cs_tkn,
              bdfe=(dfe_6pt,),
              # cdfe=(dfe_6pt,),
              )
a.fixed.cth = None
a.save()

#%% NEW BMAP PLOT
sdir = root_dir + '/result/final_files/'

suf = 'popdiv'
bv = np.load(sdir + suf + '.bvals.npy')
pi = np.load(sdir + suf + '.pivals.npy')
nn = np.load(sdir + suf + '.nvals.npy')
step = len(bv) / 100
si = np.argsort(bv)
bv, pi, nn = bv[si], pi[si], nn[si]
bmean, pimean = [], []
for i in range(0, bv.size, step):
    j = i + step
    bmean.append(np.average(bv[i:j], weights=nn[i:j]))
    pimean.append(np.average(pi[i:j], weights=nn[i:j]))
print len(pimean)
plt.plot(bmean, pimean)
plt.show()

#%% RESCALED B MAPS
bdir = root_dir + '/precalc/ape_cons95_euarchontoglires_neut30_filtered/'
fmt = 'AA_Map_10kb.{}.ape_cons95_euarchontoglires_neut30_filtered.t{}.merge5e07.bkgd'
exps = [-4.5, -4, -3.5, -3, -2.5, -2]
# plt.figure(figsize=(10, 5))
for x in exps:
    # concatenate all chroms for each t value
    # bval, size = [], []
    tv = '{:.8f}'.format(10**x)
    # for ch in human_autosomes:
    #     f = bdir + fmt.format(ch, tv)
    #     b, s = np.loadtxt(f).T
    #     bval.append(b)
    #     size.append(s)
    # bval = np.concatenate(bval)
    # size = np.concatenate(size)
    # # adjust bvals
    # bval *= np.log1p(-1.0 / 100)
    # # reset u del to 1e-08, convert to B
    # bval = np.exp(bval / 7.4)
    fsave = bdir + 'composite.{}.npy'.format(tv)
    # np.save(fsave, np.column_stack((bval, size)))
    bval, size = np.load(fsave).T
    borig = np.exp(np.log(bval) * 7.4)
    # add values to histogram
    lbl = r'$10^{%.1f}$' %x
    bi_1 = np.argsort(bval)
    bi_2 = np.argsort(borig)
    sz_1 = np.cumsum(size[bi_1])
    sz_2 = np.cumsum(size[bi_2])
    plt.step(borig[bi_2], sz_2 / sz_2[-1], label=lbl + ' u=7.4e-08')
    plt.step(bval[bi_1], sz_1 / sz_1[-1], label=lbl + ' u=1.0e-08')

    # plt.hist(bval, bins=100, label=lbl + ' u=1e-08', weights=size,
    #          histtype='step', alpha=0.9, lw=1)
    # plt.hist(borig, bins=100, label=lbl + ' u=7.4e-08', weights=size,
    #          histtype='step', alpha=1, lw=1)
    plt.legend()
    plt.savefig(fsave.replace('.npy', '.png'), dpi=256)
    plt.close()
    # plt.show()

#%% CONVERTING MCVICKER MAPS TO MATCH FORMAT
ifmt = root_dir + '/mcvicker/hg19_bkgd/{}.hg19.bkgd'
ofmt = root_dir + '/precalc/mcvicker/' \
                  'AA_Map_10kb.{}.mcvicker.t0.00100000.merge5e07.bkgd'
for ch in human_autosomes:
    b_1, seg = np.loadtxt(ifmt.format(ch), dtype=int).T
    b_2 = (np.log(np.maximum(b_1 / 1000.0, 0.001)) / np.log1p(-0.01)).astype(int)
    np.savetxt(ofmt.format(ch), np.column_stack((b_2, seg)), fmt='%d %d')

#%%
def pop_div_analysis(rdir, titl):
    fdir = root_dir + '/result/final_files/{}/'.format(rdir)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    slist = [f for f in os.listdir(fdir) if ('predsort' in f) and
             (f.endswith('.txt'))]

    f_init = fdir + flist[0]
    f_sort = fdir + slist[0]
    f_save = fdir[:-1] + '_pop_div_analysis.png'
    rst = RunStruct(init=f_init)

    plt.figure(figsize=(5, 5))
    plt.subplots_adjust(left=0.16, bottom=0.1, right=0.98, hspace=0.3,
                        top=0.90, wspace=0.35)
    tit = titl
    plt.title(tit, fontsize=16)

    if titl == 'YRI/CEU div':
        suf = 'popdiv'
        f_save = fdir[:-1] + '_pop_div_analysis.png'

    else:
        suf = 'pi'
        f_save = fdir[:-1] + '_pi_analysis.png'

    bv = np.load(fdir + suf + '.bvals.npy')
    pi = np.load(fdir + suf + '.pivals.npy')
    nn = np.load(fdir + suf + '.nvals.npy')
    step = len(bv) / 100
    si = np.argsort(bv)
    bv, pi, nn = bv[si], pi[si], nn[si]
    bmean, pimean = [], []
    for i in range(0, bv.size, step):
        j = i + step
        bmean.append(np.average(bv[i:j], weights=nn[i:j]))
        pimean.append(np.average(pi[i:j], weights=nn[i:j]))
    pimean = np.array(pimean)[2:-5]
    bmean = np.array(bmean)[2:-5]
    # print len(pimean)
    # pi0 = rst.stat.meanpi * rst.fixed.tau_init / rst.params[-1]
    # div, pi, pred = np.loadtxt(f_sort)[2:-5].T
    # b = pred / pi0
    # xi = np.arange(div.size)
    # pi /= pi0
    # if rdir == 'ceu_ape95':
    #     yri_init = root_dir +'/result/final_files/ape95_minbs_01/YRI.ape95_euarchontoglires35_filtered.BS1.6.CS0.0.iprm_02.bth_0.000.191024162056.final.composite.txt'
    #     yri_rst = RunStruct(init=yri_init)
    #     yri_pi0 = yri_rst.stat.meanpi * yri_rst.fixed.tau_init / yri_rst.params[-1]
    #     pred /= yri_pi0
    # else:
    #     pred /= pi0
    # if rdir == 'yri_ceu_pop_div':
    #     pi *= (pred.mean() / pi.mean())

    # run linear regression on pi/pred
    tslope, tintercept, lo_slope, hi_slope = theilslopes(bmean, pimean, 0.95)
    slope_ci = hi_slope-lo_slope
    slope, intercept, rvalue, pvalue, stderr = linregress(bmean, pimean)
    lmodel = r'$y = %.3ex (\pm %.3e) + %.3e$' %(slope, slope_ci, intercept, )
    lfunc = lambda x: slope*x + intercept

    # # SORTED MAP 1
    # plt.subplot(121)
    # plt.axhline(y=1, color='k', ls='--', alpha=0.8)
    # plt.text(50, 1.03, 'w/o background selection', ha='center', va='center')
    # plt.axhline(y=meanpi, color='darkslategray', ls=':', alpha=0.8)
    # plt.text(65, meanpi - 0.03, 'mean '+r'$\pi/\pi_0$', ha='left',
    #          va='center', color='darkslategray')
    # plt.plot(xi, pi, label='observed', color='darkslategray')
    # plt.plot(xi, pred, label='predicted', color='fuchsia')
    # plt.ylabel(r'$\pi/\pi_0$', fontsize=14)
    # plt.ylim(0.48, 1.2)
    # if rdir == 'mcvicker':
    #     plt.ylim(0, 1.2)
    # plt.xlabel('background selection bin')
    # plt.legend(loc='lower right')

    # SORTED MAP 2
    # plt.subplot(122)
    # plt.axhline(y=1, color='k', ls='--', alpha=0.8)
    # plt.text(50, 1.03, 'w/o background selection', ha='center', va='center')
    pimin = pimean.min()
    pimax = pimean.max()
    plt.text(0.8, 1.1*pimax, lmodel, ha='center', va='center', color='k')
    # plt.axhline(y=meanpi, color='darkslategray', ls=':', alpha=0.8)
    plt.plot(bmean, pimean, label='observed', color='darkslategray')
    # plt.plot(b, pred, label='predicted', color='fuchsia')
    # plt.plot([0, 1], [0, 1], label=r'$y=x$', ls='--', lw=1, color='b')
    plt.plot([0.5, 1.05], map(lfunc, [0.5, 1.05]), label='linear fit',
             ls='--', lw=1, color='b')
    plt.ylabel(r'$\pi/\pi_0$', fontsize=14)
    plt.ylabel('scaled pi')
    plt.ylim(0.9*pimin, 1.2*pimax)
    plt.xlabel('B value')
    plt.xlim(0.5, 1.05)
    if rdir == 'mcvicker':
        plt.ylim(0, 1.2)
        plt.xlim(0, 1.05)
    plt.legend(loc='lower right')

    plt.savefig(f_save, dpi=256)
    plt.close()

