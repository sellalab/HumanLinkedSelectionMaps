# -*- coding: utf-8 -*-

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:42:13 2018

@author: davidmurphy
"""

from classes.bkgdmap import BkgdMapReader
from classes.runstruct import ChromStruct
#import numpy as np
import os

#%%
bbf = '/Users/davidmurphy/GoogleDrive/linked_selection/precalc/std_split/' \
    'AA_Map_10kb.chr14.primate_cons95_Segments.t0.00010000.merge1e07.bkgd'
smode = True
br = BkgdMapReader(bbf, safemode=smode)
bm = br.get_bmap(safemode=smode)

#%%
#def merge_bmaps(cst, coef, plen):
# TODO: move partition files and compress them
# collect ordered b values and segment lengths

ch = 'chr2'
coef = 10**-4.5
plen = 1e7
anno = 'primate_cons95_Segments'
cst = ChromStruct(chrom=ch, bdir='std_split_pts',bscl=1e2,
                  tkn='pr95.cleanrun')

# create the new suffix and new directory for finished partitions    
n = int(cst.chlen / plen) + 1 # num partitions of chrom
#plen_str = '{:.0e}'.format(plen).replace('+', '')

# create new file for finished merge
merge_bfile = cst.bkgd_file(anno, coef, plen=plen, merge=1)
# make a new directory to put partition files into
#partition_dir = '{}/partitions{}'.format(cst.bmap_dir, plen_str)
#if not os.path.isdir(partition_dir):
#    os.mkdir(partition_dir)
#x = os.popen('tail -n 10 {}'.format(merge_bfile)).read()

# get the list of partition files for the current chrom
#f_list = [cst.bkgd_file(anno, coef, pidx, plen) for pidx in xrange(n)]

# =============================================================================
# 
# import gzip
# import shutil
# with open('file.txt', 'rb') as f_in, gzip.open('file.txt.gz', 'wb') as f_out:
#     shutil.copyfileobj(f_in, f_out)
# =============================================================================

#%%
# create a reader to read each of the files
br = BkgdMapReader()
br.merge_partitions(f_list)
br.write_bfile()

#%%
bmp = br.get_bmap()
print 'ok;'
#%%
# move partition files into their new subdirectory
#for f in f_list:
#    f_new = f.replace(cst.bmap_dir, partition_dir)
##    os.rename(f, f_new)

#%%
# get the merged map and write it to a new file if there are no errors
#bm = br.get_bmap()
#
##%%
#if bm.segs.sum() != cst.chlen:
#    # unified maps must sum to exact chrom length    
#    print 'WARNING: {}={}bp, partition={}bp'.format(
#            cst.chrom, cst.chlen, bm.segs.sum())
#else:
#    bscale = br.info['BKGD_SCALE']
#    bm.write_bmap(merge_bfile, bscale, info=br.info)
#        for f in ls:
#            os.
#return bm, br

#%%


#bm, br = merge_bmaps(cst, coef, plen)
#n = int(cst.chlen / plen) + 1 # num partitions of chrom
#br = BkgdMapReader()
#
#ls = [cst.bkgd_file(anno, t, pidx, plen) for pidx in xrange(n)]
#br.merge_bfiles(ls)
#bm = br.get_bmap()

#s0 = ls[0].replace('p0.1e07', '1e07.merged')
#bm.write_bmap(s0, 100)
       
#bm, br = merge_bmaps(cst, t, 1e6)

#%%

# =============================================================================
# for ch in cst.chroms[7:10]:  
#     cst.chrom = ch    
#     bmg, brg = merge_bmaps(cst, t, 1e6)
#     
#     #info = br.info
#     #fout = '/Users/davidmurphy/Desktop/test.bkgd'
#     #bm.write_bmap(fout, sc, info)
#     
#     br = BkgdMapReader(cst.bkgd_file('primate_cons95_Segments', t))
#     bm = br.get_bmap()
#     
#     # load random points and compare both sets of b values
#     bx = np.load(cst.rand_file(25000))['randneut']
#     bm_x, bmg_x = bm.interp_bkgd(bx), bmg.interp_bkgd(bx)
#     
#     ba = np.column_stack((bx, bm_x, bmg_x))    
# 
#     bl.append(ba)
# =============================================================================

#%%
import matplotlib
font = {'size': 24}
matplotlib.rc('font', **font)
plt = matplotlib.pyplot
import seaborn

#%%
plt.figure()
plt.hist(diff, bins=500)
plt.figure()
plt.hist(diff, bins=50, cumulative=1, histtype='step', normed=True, lw=2)
#plt.figure(figsize=(10,10))
#plt.scatter(bx, diff, color='orange')
#plt.figure(figsize=(12, 6))
#plt.plot(bx, bmg_x, label='merged', lw=1, color='purple')
#plt.plot(bx, bm_x, label='single', lw=1, color='darkorange', alpha=0.8)
#plt.legend()
#plt.show()np.array([0])
#np.array([0])
