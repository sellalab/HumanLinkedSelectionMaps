__author__ = 'davidmurphy'

import seaborn
import numpy as np
from itertools import izip
import matplotlib.pyplot as plt
from classes.bkgdmap import BkgdMapReader
from collections import defaultdict, Counter
from classes.runstruct import ChromStruct, root_dir, izip


#%%
bdir = 'std_split_pts'
cst = ChromStruct(chrom='chr1', bdir=bdir)
an = cst.bs_annos[0]
pl = 1e7

#%%


#%%
datadict = defaultdict(Counter)
for ti in cst.bdfe[0]:
# for ti in [1e-02]:
    for ch in cst.chroms:
        cst.chrom = ch
        bf = cst.bkgd_file(an, ti, plen=pl, merge=True)
        bmr = BkgdMapReader(bf)
        bmp = bmr.get_bmap()
        for s, b in izip(bmp.segs, bmp.bkgd):
            datadict[ti][b] += s

#%%
lbl = [r'$10^{%.1f}$' % x for x in np.linspace(-4.5, -2.0, 6)]
for ti in cst.bdfe[0]:
    cntr = datadict[ti]
    tot = sum(cntr.values())
    bv = np.array(sorted(cntr.keys()))
    sg = np.cumsum([cntr[k] / tot for k in bv])
    plt.step(bv / -7.4e-08, sg, label=str(ti), lw=0.75)
fsave = '/Users/davidmurphy/Desktop/logcumplot.png'
plt.xscale('log')
plt.ylabel('cumulative fraction of b values')
plt.xlabel('-b values (udel=1)')
plt.legend(lbl)
plt.savefig(fsave, dpi=256)
# plt.show()

#%%