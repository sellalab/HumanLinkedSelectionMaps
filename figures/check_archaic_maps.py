__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from classes.runstruct import root_dir, cst_from_fldr, chromosome_length, \
    human_autosomes, ChromStruct
from classes.phylotree import parse_exptotsub
from data_processing.data_tools import cg_hypermutable
from figures.common_functions import format_panels, get_bbins, \
    get_telomere_dict, get_centromere_dict, predict_loess

# u0 for scaling parameter estimates
u0 = 1.4e-08
final_dir = root_dir + '/result/final_files'

#%% GET DICTIONARY OF NEUTRAL MASKS FOR EACH CHROM
def get_nmsk(chrom):
    nstr = root_dir + '/data/nmsk/{}.euarchontoglires.0.35.gmsk0.1.nmsk.npz'
    nmsk = np.load(nstr.format(ch))['neutmask']

    return nmsk

ndict = {}
for ch in human_autosomes:
    nm = get_nmsk(ch)
    ndict[ch] = nm


#%% FUNCTION FOR LOADING BMAP AS AN ARRAY OF VALUES FOR EVERY POSITION IN CHROM
def get_bmap_array(fldr, chrom):
    # set bmap file path and load bmap
    f_bmap = final_dir + '/{}/bmap/{}.bmap.txt'.format(fldr, chrom)
    bmap = np.loadtxt(f_bmap, dtype=int)

    # create empty array to fill with bmap values
    barr = np.zeros(shape=chromosome_length(chrom), dtype='u2')

    # fill entire chrom length array with bmap values
    i = 0
    for (b, l) in bmap:
        barr[i:i+l] = b
        i += l

    # check that the full chrom was covered
    assert i == barr.size

    return barr


fldr = 'cadd94_gmask_mnb_378'
bdict = {}
for ch in human_autosomes:
    ba = get_bmap_array(fldr, ch)
    bdict[ch] = ba


#%% FUNCTION THAT LOADS COORDS FROM AFR ARCHAIC MAPS FOR EACH CHROM TO A DICT
def get_afr_archaic(pop):
    # set archie calls file path
    fdir = root_dir + '/data/snps/ArchIE-calls-YRI-MSL/'
    f_archaic = fdir + '{}-freq-map.bed'.format(pop)

    # load archie calls and store rows in dict using chroms as keys
    adict = defaultdict(list)
    with open(f_archaic, 'r') as f:
        for line in f:
            ch, start, end, prob = line.split()
            adict[ch].append((int(start), int(end), float(prob)))

    return adict


adict = get_afr_archaic('YRI')


#%% FUNCTION THAT LOADS COORDS FROM EUR ARCHAIC MAPS FOR EACH CHROM TO A DICT
def get_eur_archaic():
    # set archie calls file path
    f_archaic = root_dir + '/data/snps/archaic.segments.txt'
    # load archie calls and store rows in dict using chroms as keys
    edict = defaultdict(list)
    with open(f_archaic, 'r') as f:
        f.next()
        for line in f:
            line = line.split()
            if line[-1] != 'HMM':
                continue
            ch, start, end = line[1:4]
            ch = 'chr' + ch
            prob = line[9]
            if line[8] == 'SouthAsia':
            # if line[8] == 'EastAsia':
            # if line[8] == 'WestEurasia':
            # if line[7] == 'Abkhazia':
                edict[ch].append((int(start), int(end), float(prob)))

    return edict


edict = get_eur_archaic()


#%% PAIR EACH SEGMENT CALLED FOR ARCHAIC WITH ITS MEAN B VALUE
def archaic_bmap_paired(ndict, bdict, adict):
    probs, bvals = [], []
    for ch in human_autosomes:
        ba = bdict[ch]
        nm = ndict[ch]
        for (start, end, prob) in adict[ch]:
            b_block = ba[start:end]
            n_block = nm[start:end]
            mean_b = np.mean(b_block[n_block])
            probs.append(prob)
            bvals.append(mean_b)

    return probs, bvals


#%%
# pr, bv = archaic_bmap_paired(ndict, bdict, edict)
pr, bv = archaic_bmap_paired(ndict, bdict, adict)

#
pr = np.array(pr)
bv = np.array(bv)
pr = pr[~np.isnan(bv)]
bv = bv[~np.isnan(bv)]
#
sidx = np.argsort(bv)
bv = bv[sidx]
pr = pr[sidx]
xt = np.arange(min(bv), max(bv), 1)
plo = predict_loess(bv, pr, None, 0.1, xt)


#
mean_b, mean_p = [], []
step = len(bv) / 1000
for i in range(0, len(bv), step):
    mean_b.append(np.mean(bv[i:i+step]))
    mean_p.append(np.mean(pr[i:i+step]))


#
plt.figure()
plt.scatter(mean_b, mean_p, alpha=0.5, s=5, color='deepskyblue')
plt.plot(xt, plo, color='cornflowerblue')
# f_save = final_dir + '/sfigs/b_vs_EUR_archaic.png'
# f_save = final_dir + '/sfigs/b_vs_PAK_archaic.png'
f_save = final_dir + '/sfigs/b_vs_AFR_archaic.png'

plt.savefig(f_save, dpi=512)
plt.close()


#%%
plt.figure()
plt.hist(bv, bins=100)
f_save = final_dir + '/sfigs/b_vs_EUR_archaic_hist.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%%
plt.figure()
m = (pr > 0)
xt = np.arange(min(bv), max(bv), 1)
plo = predict_loess(bv[m], pr[m], None, 0.1, xt)
# plt.scatter(bv[m], pr[m], alpha=0.5, s=6, color='dodgerblue')
plt.plot(xt, plo, color='blue')
f_save = final_dir + '/sfigs/b_vs_archaic_mask.png'
plt.savefig(f_save, dpi=512)
plt.close()

#%%