__author__ = 'davidmurphy'


import os
import numpy as np
from datetime import datetime as dt
from sys import argv, stdout, stderr
from classes.geneticmap import GeneticMap
from classes.annopoints import AnnoPoints
from classes.cscalculator import CSCalculator
from data_processing.data_tools import mask_segments
from classes.runstruct import ChromStruct, root_dir, chromosome_length


def annotated_hc_substitutions(ch, anno):
    """get positions of HC-substitutions in specific anno"""
    # set cst for filepaths
    cst = ChromStruct(ch)

    # get filtered cons mask
    f_con = cst.bs_target(anno)
    c_seg = np.loadtxt(f_con, usecols=(1, 2), dtype=int)
    msk = np.zeros(cst.chlen, dtype=bool)
    msk = mask_segments(msk, c_seg, flipoff=False)

    # get indices of HC derived substitutions
    f_der = cst.data + '/div/hcder/{}.hc.anc.der.idx.npz'.format(ch)
    d_pos = np.load(f_der)['derived'][:,0].astype(int)

    # get derived positions in the annotations (d_pos where msk is true)
    c_pos = d_pos[msk[d_pos]]

    # print out the number of total and annotated derived sites
    msg = 'total/annotated {} {}'.format(len(d_pos), len(c_pos))
    stderr.write(msg)
    stdout.flush()

    c_an = '{}_hc_subs'.format(anno)
    fout = cst.cs_target(c_an)
    pth = '/'.join(fout.split('/')[:-1])
    if not os.path.isdir(pth):
        os.mkdir(pth)
    np.savez_compressed(fout, pos=c_pos)

    return None


def exnex_conserved_substitutions(chrom):
    # set annotations
    x_an = 'ape_cons94_exonic'
    nx_an = 'ape_cons94_nonexonic'

    # set paths to files
    f_path = root_dir + '/data/bsanno'
    f_ex = f_path + '/{a}/{c}.{a}.bed'.format(a=x_an, c=chrom)
    f_nex = f_path + '/{a}/{c}.{a}.bed'.format(a=nx_an, c=chrom)
    f_der = root_dir + '/data/div/hcder/{}.hc.anc.der.idx.npz'.format(chrom)

    # initialize masks for ex/nex annotations
    xmsk = np.zeros(shape=chromosome_length(chrom), dtype=bool)
    nxmsk = np.zeros(shape=chromosome_length(chrom), dtype=bool)

    # group files, masks in lists
    f_in = [f_ex, f_nex]
    msks = [xmsk, nxmsk]

    for (f, m) in zip(f_in, msks):
        c_seg = np.loadtxt(f, usecols=(1, 2), dtype=int)
        mask_segments(m, c_seg, flipoff=False)

    # get indices of HC derived substitutions
    d_pos = np.load(f_der)['derived'][:, 0].astype(int)

    # get conserved derived exonic positions (d_pos where conserved is true)
    x_pos = d_pos[xmsk[d_pos]]
    nx_pos = d_pos[nxmsk[d_pos]]
    points = [x_pos, nx_pos]

    # set file save paths, output annotations, create directories if needed
    xa = 'ape_cons94_hc_der_exonic'
    nxa = 'ape_cons94_hc_der_nonexonic'
    annos = [xa, nxa]
    s_path = root_dir + '/data/csanno'
    for (an, pt) in zip(annos, points):
        pth = '{}/{}'.format(s_path, an)
        if not os.path.isdir(pth):
            os.mkdir(pth)
        f_save = pth + '/{}.{}.npz'.format(chrom, an)
        np.savez_compressed(f_save, pos=pt)

    return None


def ns_nonns_conserved_substitutions(chrom, anno):

    # set paths to nonsyn and HC-derived position files
    cs_path = root_dir + '/data/csanno'
    f_ns = cs_path + '/nonsyn/{}.nonsyn.npz'.format(chrom)
    f_hc = root_dir + '/data/div/hcder/{}.hc.anc.der.idx.npz'.format(chrom)

    # set path to conservation annotation file for anno
    bs_path = root_dir + '/data/bsanno'
    f_an = bs_path + '/{a}/{c}.{a}.bed'.format(a=anno, c=chrom)

    # load nonsyn and HC sites, make sure nonsyn is a complete subset of HC
    ns = np.load(f_ns)['pos'].astype(int)
    hc = np.load(f_hc)['derived'][:, 0].astype(int)
    nsidx = np.in1d(hc, ns)
    assert np.sum(nsidx) == ns.size

    # create a subset of NOT nonsyn HC sites
    nonns = hc[~nsidx]

    # create a mask of annotated sites
    mask = np.zeros(shape=chromosome_length(chrom), dtype=bool)
    an = np.loadtxt(f_an, usecols=(1, 2), dtype=int)
    mask_segments(mask, an, flipoff=False)

    # get positions of conserved NS and conserved non-NS sites
    ns_pos = ns[mask[ns]]
    nonns_pos = nonns[mask[nonns]]
    points = [ns_pos, nonns_pos]

    # set file save paths, output annotations, create directories if needed
    ns_an = 'nonsyn_{}_hc_subs'.format(anno)
    nonns_an = 'other_{}_hc_subs'.format(anno)
    annos = [ns_an, nonns_an]
    for (an, pt) in zip(annos, points):
        pth = '{}/{}'.format(cs_path, an)
        if not os.path.isdir(pth):
            os.mkdir(pth)
        f_save = pth + '/{}.{}.npz'.format(chrom, an)
        np.savez_compressed(f_save, pos=pt)

    return None


def ns_nonns_conserved_substitutions_sampled(chrom, anno, pop, sample):

    # set paths to nonsyn and HC-derived position files
    cs_path = root_dir + '/data/csanno'
    f_ns_str = '/{p}_nonsyn_s{s}/{c}.{p}_nonsyn_s{s}.npz'
    f_ns = cs_path + f_ns_str.format(p=pop, s=sample, c=chrom)
    f_hc_str = '/data/div/hcder/{}.{}.HC.subst.sample.{}.npz'
    f_hc = root_dir + f_hc_str.format(chrom, pop, sample)

    # set path to conservation annotation file for anno
    bs_path = root_dir + '/data/bsanno'
    f_an = bs_path + '/{a}/{c}.{a}.bed'.format(a=anno, c=chrom)

    # load nonsyn and HC sites, make sure nonsyn is a complete subset of HC
    ns = np.load(f_ns)['pos'].astype(int)
    hc = np.load(f_hc)['sub'][:, 0].astype(int)
    nsidx = np.in1d(hc, ns)
    assert np.sum(nsidx) == ns.size

    # create a subset of NOT nonsyn HC sites
    nonns = hc[~nsidx]

    # create a mask of annotated sites
    mask = np.zeros(shape=chromosome_length(chrom), dtype=bool)
    an = np.loadtxt(f_an, usecols=(1, 2), dtype=int)
    mask_segments(mask, an, flipoff=False)

    # get positions of conserved NS and conserved non-NS sites
    ns_pos = ns[mask[ns]]
    nonns_pos = nonns[mask[nonns]]
    points = [ns_pos, nonns_pos]

    # print message about the number of sites in each category
    msg = 'NS={} OTHER={}\n'.format(len(ns_pos), len(nonns_pos))
    stderr.write(msg)
    stdout.flush()

    # set file save paths, output annotations, create directories if needed
    ns_an = 'nonsyn_{}_{}_s{}_hc_subs'.format(anno, pop, sample)
    nonns_an = 'other_{}_{}_s{}_hc_subs'.format(anno, pop, sample)
    annos = [ns_an, nonns_an]
    for (an, pt) in zip(annos, points):
        pth = '{}/{}'.format(cs_path, an)
        if not os.path.isdir(pth):
            os.mkdir(pth)
        f_save = pth + '/{}.{}.npz'.format(chrom, an)
        np.savez_compressed(f_save, pos=pt)

    return None


# ns_nonns_conserved_substitutions('chr22', 'cadd94')


def build_csmaps(f_init, ch, expo):
    """build sweeps map based for given chrom and s value exponent"""
    cst = ChromStruct(chrom=ch, init=f_init)
    coef = 10 ** expo
    ne = 2e4

    # run exonic and nonexonic in same job
    for an in cst.cs_annos:
        # file paths and params
        # cst = ChromStruct(ch, fcl=an, cdir=cdir)
        cfile = cst.cmap_file(an, coef)
        afile = cst.cs_target(anno=an)

        # annotation and map classes
        gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)
        ano = AnnoPoints(gmp, afile)

        # initialize CS calculator and save grid of fixation rates
        csc = CSCalculator(ano, coef, ne)

        start = dt.now()
        # create path to precalc files if it doesnt exist
        c_pth = '/'.join(cfile.split('/')[:-1])
        if not os.path.isdir(c_pth):
            os.mkdir(c_pth)
        c = csc.grid_fixrates(lim=True, fout=cfile)[0]
        t = dt.now() - start
        msg = 'cmap length: {} segments\nruntime: {}\n'.format(len(c), t)
        stderr.write(msg)
        stdout.flush()


def main():
    if len(argv) == 3:
        ch = argv[1]
        an = argv[2]
        # annotated_hc_substitutions(ch, an)
        ns_nonns_conserved_substitutions(ch, an)
    elif len(argv) == 4:
        f_init = argv[1]
        ch = argv[2]
        expo = -float(argv[3])
        build_csmaps(f_init, ch, expo)
    elif len(argv) == 5:
        ch, an, pop, ss = argv[1:]
        ns_nonns_conserved_substitutions_sampled(ch, an, pop, ss)
    else:
        print 'usage_1: cs_annotations <chrom> <anno>'
        print 'usage_2: cs_annotations <init> <chrom> <expo>'
        print 'usage_3: cs_annotations <chrom> <anno> <pop> <sample>'
        exit()


if __name__ == '__main__':
    main()
