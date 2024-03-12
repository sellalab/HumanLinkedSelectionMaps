__author__ = 'davidmurphy'


import os
import seaborn
import numpy as np
from sys import argv
from itertools import izip
import matplotlib.pyplot as plt
from classes.runstruct import ChromStruct, root_dir
from precalc.lh_inputs import load_saved, adjust_arrays

class NeData:
    def __init__(self, txt_file):
        # create list for raw values
        low, high, val = [], [], []
        with open(txt_file, 'r') as f:
            for line in f:
                lo, hi, v = line.split('\t')
                low.append(float(lo))
                high.append(float(hi))
                val.append(float(v))
        # create arrays from value lists
        self.low = np.array(low)
        self.high = np.array(high)
        self.val = np.array(val)

    @property
    def dtime(self):
        return self.high - self.low

    def tmrca(self, b):
        t = self.dtime / 30.0
        n = self.val
        jr = range(len(t) - 1)

        exp = np.exp(-(t / (2*n)) * b)
        prod = np.product([exp[j] for j in jr])
        return np.sum(prod * (1 - exp) * 2 * n * b)


def binned_b_and_tmrca(cum_pos, ms, bs, scale, slide=None):
    # get TMRCA for each B value
    # f = '/Users/davidmurphy/Desktop/Ne_T.txt'
    f = '/ifs/data/c2b2/gs_lab/dam2214/run/Ne_T.txt'
    ne = NeData(f)
    tmrca = np.array([ne.tmrca(b) for b in bs])

    # create indices for sliding windows
    if slide:
        assert slide < 1
        windows = np.arange(0, cum_pos[-1], scale * slide)
        si = np.searchsorted(cum_pos[ms], windows)
        idx = si[:-2]
        jdx = si[2:]
    # create indices for nonoverlapping windows
    else:
        # create upper, lower indices to sort data into genomic windows
        windows = np.arange(0, cum_pos[-1], scale)
        jdx = np.searchsorted(cum_pos[ms], windows)
        idx = np.concatenate(([0], jdx[:-1]))

    # prepare empty lists for each relevant statisic
    bv, tm = [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            bv.append(bs[i:j].mean())
            tm.append(tmrca[i:j].mean())

    bv = np.array(bv)
    tm = np.array(tm)

    return np.column_stack((bv, tm))


def do_binning():
    if len(argv) != 2:
        print 'usage: calc_rsq <folder_name>'
        exit(1)
    scale = 1e6
    slide = 0.5
    chrom = 'chr1'

    fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst, chroms=[chrom])

    # mask and rescale maps
    bs, cs, nu, nt, dv, pl, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)

    # convert segments into positions for r r squared calc
    cum_pos = np.cumsum(sg)

    params = cst.params
    uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
    bwt = uvec / cst.fixed.u_fix
    bsx = np.exp(np.dot(bs, bwt))

    # get predicted and observed diversity over the chromosome
    po_arr = binned_b_and_tmrca(cum_pos, msk, bsx, scale, slide)

    # save predicted and observed diversity to file
    fpth = '/'.join(f_init.split('/')[:-1])
    if slide:
        fmt = fpth + '/{}.{:.2e}win_{:.1f}slide_B_and_TMRCA.txt'
        fout = fmt.format(chrom, scale, slide)
    else:
        fmt = fpth + '/{}.{:.2e}win_B_and_TMRCA.txt'
        fout = fmt.format(chrom, scale)

    np.savetxt(fout, po_arr)


def do_plotting():
    f = root_dir + '/result/final_files/ape_cons95_clean/chr1.1.00e+06win_0.5slide_B_and_TMRCA.txt'
    b, t = np.loadtxt(f).T
    xi = np.arange(0, b.size / 2.0, 0.5)
    plt.plot(xi, b / b.mean(), label='B')
    plt.plot(xi, t / t.mean(), label='E[TMRCA]')
    plt.xticks(range(25, 250, 50))
    plt.xlabel('chr1 position (Mb)')
    plt.yticks(np.arange(0.7,1.21,0.05))
    plt.ylabel(r'$y/\bar{y}$')
    plt.ylim(0.67, 1.22)
    plt.legend()
    f_fig = f.replace('.txt', '.png')
    plt.savefig(f_fig, dpi=256)
    plt.close()
    # plt.show()


def main():
    # do_binning()
    do_plotting()


if __name__ == '__main__':
    main()
