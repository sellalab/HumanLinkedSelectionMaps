__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from classes.runstruct import ChromStruct, RunStruct
from classes.geneticmap import GeneticMap

spec = 'ape primate prosimian euarchontoglires laurasiatheria ' \
       'mammal fish'.split()


def genmap_bin_cons(ch, bin_size, flag='depth'):
    """count conserved sites in genetic map unit bins on different scales"""
    # init ChromStruct and GeneticMap classes
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)

    # create bins for counting up conserved sites in each segment
    gmin, gmax = gmp.gmap_limits
    bins = np.arange(gmin, gmax + bin_size, bin_size)

    # get conserved segments from each species and store genetic positions
    cnt = []
    pcts = map(str, range(91, 99))
    if flag == 'depth':
        cols = spec
        for spc in spec:
            # load segments for annotation
            anno = '{}_cons94_gmask'.format(spc)
            f_cons = cst.bs_target(anno)
            # convert segments into array of genetic positions
            pos = []
            for (start, end) in np.loadtxt(f_cons, usecols=(1, 2)):
                pos.append(np.arange(start, end, 1))
            pos = np.concatenate(pos)
            gpos = gmp.interp_gpos(pos)
            # sort segment positions into bins
            cnt.append(np.histogram(gpos, bins=bins)[0])
    else:
        assert flag == 'percent'
        cols = pcts
        for pc in pcts:
            # load segments for annotation
            anno = 'fish_cons{}_gmask'.format(pc)
            f_cons = cst.bs_target(anno)
            # convert segments into array of genetic positions
            pos = []
            for (start, end) in np.loadtxt(f_cons, usecols=(1, 2)):
                pos.append(np.arange(start, end, 1))
            pos = np.concatenate(pos)
            gpos = gmp.interp_gpos(pos)
            # sort segment positions into bins
            cnt.append(np.histogram(gpos, bins=bins)[0])

    # save the counts per bin for each species
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_save = '{}/{}.{}.gmap.dist.bin_{:.2e}.txt'.format(pth, ch, flag, bin_size)
    # l_fmt = '{} {:.0f} {:.0f} {:.0f}\n'
    l_fmt = '{} {}\n'
    with open(f_save, 'w') as f:
        f.write('#bin ' + ' '.join(cols) + '\n')
        for i in range(len(bins)-1):
            # l_a, l_e, l_f = cnt[0][i], cnt[1][i], cnt[2][i]
            data = [c[i] for c in cnt]
            if all(c == 0 for c in (data)):
                continue
            else:
                d_line = ' '.join('{:.0f}'.format(c) for c in data)
                line = l_fmt.format(bins[i], d_line)
                f.write(line)

    return None


def genmap_compare_cadd(ch, bin_size, flag='depth'):
    """count conserved sites in genetic map unit bins on different scales"""
    # init ChromStruct and GeneticMap classes
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)

    # create bins for counting up conserved sites in each segment
    gmin, gmax = gmp.gmap_limits
    bins = np.arange(gmin, gmax + bin_size, bin_size)

    # get conserved segments from each species and store genetic positions
    cnt = []
    # f_cons = cst.bs_target('ape_cons94_clean')
    # f_cadd = cst.bs_target('cadd94')
    # f_cons = cst.bs_target('fish_cons94_gmask')
    f_cons = cst.bs_target('ape_cons94_gmask')
    f_cadd = cst.bs_target('cadd94_gmask')
    for f_in in (f_cons, f_cadd):
        # convert segments into array of genetic positions
        pos = []
        for (start, end) in np.loadtxt(f_in, usecols=(1, 2)):
            pos.append(np.arange(start, end, 1))
        pos = np.concatenate(pos)
        gpos = gmp.interp_gpos(pos)
        # sort segment positions into bins
        cnt.append(np.histogram(gpos, bins=bins)[0])

    # save the counts per bin for each species
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_save = '{}/{}.{}.gmap.dist.bin_{:.2e}.txt'.format(pth, ch, flag, bin_size)
    l_fmt = '{} {}\n'
    cols = 'ape cadd'
    with open(f_save, 'w') as f:
        f.write('#bin ' + ' '.join(cols) + '\n')
        for i in range(len(bins)-1):
            data = [c[i] for c in cnt]
            if all(c == 0 for c in (data)):
                continue
            else:
                d_line = ' '.join('{:.0f}'.format(c) for c in data)
                line = l_fmt.format(bins[i], d_line)
                f.write(line)

    return None


def genmap_compare_exnex(ch, bin_size, flag='depth'):
    """count conserved sites in genetic map unit bins on different scales"""
    # init ChromStruct and GeneticMap classes
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)

    # create bins for counting up conserved sites in each segment
    gmin, gmax = gmp.gmap_limits
    bins = np.arange(gmin, gmax + bin_size, bin_size)

    # get conserved segments from each species and store genetic positions
    cnt = []
    # f_cn = cst.bs_target('ape_cons94_clean')
    # f_ex = cst.bs_target('ape_cons94_exonic')
    # f_nex = cst.bs_target('ape_cons94_nonexonic')
    f_cn = cst.bs_target('fish_cons94_gmask')
    f_ex = cst.bs_target('fish_cons94_gmask_exonic')
    f_nex = cst.bs_target('fish_cons94_gmask_nonexonic')
    for f_in in (f_cn, f_ex, f_nex):
        # convert segments into array of genetic positions
        pos = []
        for (start, end) in np.loadtxt(f_in, usecols=(1, 2)):
            pos.append(np.arange(start, end, 1))
        pos = np.concatenate(pos)
        gpos = gmp.interp_gpos(pos)
        # sort segment positions into bins
        cnt.append(np.histogram(gpos, bins=bins)[0])

    # save the counts per bin for each species
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_save = '{}/{}.{}.gmap.dist.bin_{:.2e}.txt'.format(pth, ch, flag, bin_size)
    l_fmt = '{} {}\n'
    cols = 'all exonic other'.split()
    with open(f_save, 'w') as f:
        f.write('#bin ' + ' '.join(cols) + '\n')
        for i in range(len(bins)-1):
            data = [c[i] for c in cnt]
            if all(c == 0 for c in (data)):
                continue
            else:
                d_line = ' '.join('{:.0f}'.format(c) for c in data)
                line = l_fmt.format(bins[i], d_line)
                f.write(line)

    return None


def genmap_compare_bscs(ch, bin_size, flag='depth'):
    """count conserved sites in genetic map unit bins on different scales"""
    # bs_anno = 'cadd93'
    bs_anno = 'fish_cons94_gmask'
    # cs_anno = bs_anno + '_hc_subs'

    # init ChromStruct and GeneticMap classes
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)

    # create bins for counting up conserved sites in each segment
    gmin, gmax = gmp.gmap_limits
    bins = np.arange(gmin, gmax + bin_size, bin_size)

    # get annotation segments and store genetic positions per site
    cnt = []
    pos = []
    f_bs = cst.bs_target(bs_anno)
    for (start, end) in np.loadtxt(f_bs, usecols=(1, 2)):
        pos.append(np.arange(start, end, 1))
    pos = np.concatenate(pos)
    gpos = gmp.interp_gpos(pos)
    # sort segment positions into bins
    cnt.append(np.histogram(gpos, bins=bins)[0])

    # # get conserved HC substitution positions
    # f_cs = cst.cs_target(cs_anno)
    # pos = np.load(f_cs)['pos']
    # gpos = gmp.interp_gpos(pos)
    # cnt.append(np.histogram(gpos, bins=bins)[0])

    # get conserved HC substitution positions for each % cons
    pct = range(91, 99)
    for p in pct:
        cs_anno = 'fish_cons{}_gmask_hc_subs'.format(p)
        # cs_anno = 'ape_cons{}_clean_hc_subs'.format(p)
        f_cs = cst.cs_target(cs_anno)
        pos = np.load(f_cs)['pos']
        gpos = gmp.interp_gpos(pos)
        cnt.append(np.histogram(gpos, bins=bins)[0])

    # save the counts per bin for each species
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_save = '{}/{}.{}.gmap.dist.bin_{:.2e}.txt'.format(pth, ch, flag, bin_size)
    l_fmt = '{} {}\n'
    # cols = 'BS ' + ' '.join('{}')
    header = '#bin BS ' + ' '.join(str(p) for p in pct) + '\n'
    with open(f_save, 'w') as f:
        f.write(header)
        for i in range(len(bins)-1):
            data = [c[i] for c in cnt]
            if all(c == 0 for c in (data)):
                continue
            else:
                d_line = ' '.join('{:.0f}'.format(c) for c in data)
                line = l_fmt.format(bins[i], d_line)
                f.write(line)

    return None


def genmap_compare_cscs(ch, bin_size, flag='depth'):
    """count conserved sites in genetic map unit bins on different scales"""
    cs_anno_1 = 'cadd93_hc_subs'
    cs_anno_2 = 'ape_cons93_clean_hc_subs'

    # init ChromStruct and GeneticMap classes
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)

    # create bins for counting up conserved sites in each segment
    gmin, gmax = gmp.gmap_limits
    bins = np.arange(gmin, gmax + bin_size, bin_size)

    # get annotation segments and store genetic positions per site
    cnt = []
    for anno in [cs_anno_1, cs_anno_2]:
        f_an = cst.cs_target(anno)
        pos = np.load(f_an)['pos']
        gpos = gmp.interp_gpos(pos)
        # sort segment positions into bins
        cnt.append(np.histogram(gpos, bins=bins)[0])

    # save the counts per bin for each species
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_save = '{}/{}.{}.gmap.dist.bin_{:.2e}.txt'.format(pth, ch, flag, bin_size)
    l_fmt = '{} {}\n'
    cols = 'CADD ape'
    with open(f_save, 'w') as f:
        f.write('#bin ' + cols + '\n')
        for i in range(len(bins)-1):
            data = [c[i] for c in cnt]
            if all(c == 0 for c in (data)):
                continue
            else:
                d_line = ' '.join('{:.0f}'.format(c) for c in data)
                line = l_fmt.format(bins[i], d_line)
                f.write(line)

    return None


def main():
    if len(argv) != 4:
        print 'usage: conserved_autocorrelation <ch> <binsize> <flag>'
        exit(1)

    ch = argv[1]
    binsize = eval(argv[2])
    flag = argv[3]
    if flag == 'cadd':
        genmap_compare_cadd(ch, binsize, flag)
    elif flag == 'exnex':
        genmap_compare_exnex(ch, binsize, flag)
    elif flag == 'bscs':
        genmap_compare_bscs(ch, binsize, flag)
    elif flag == 'cscs':
        genmap_compare_cscs(ch, binsize, flag)
    else:
        genmap_bin_cons(ch, binsize, flag)


if __name__ == '__main__':
    main()