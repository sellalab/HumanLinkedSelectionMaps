__author__ = 'davidmurphy'


import os
import re
import subprocess
import numpy as np
from sys import stderr, stdout, argv
from gzip import open as zopen
from collections import Counter
from functions import relative_error
from datetime import datetime as dtime
from itertools import izip, combinations
from classes.runstruct import ChromStruct
import data_tools as dtl
from classes.phylotree import parse_tree, TreeDict


def relerr_2(cst, win):
    assert isinstance(cst, ChromStruct)
    base_count, sub_count = [], []
    for ch in cst.chroms:
        nmsk = np.load(cst.neut_masks)['neutmask']
        sbcnt = dtl.substitution_counts(nmsk, win, ch)
        bases, subs = sbcnt[:, 1:].T

        print 'boo!'


def err_msg(errstring):
    stderr.write(errstring + '\n')
    stdout.flush()


def readall(f_name):
    if f_name.endswith('.gz'):
        with zopen(f_name, 'r') as f:
            return f.read()
    else:
        with open(f_name, 'r') as f:
            return f.read()


def split_alignments(ma_fasta):
    """split a multiple alignment FASTA file and return dict of name:seq"""
    # set FASTA label regex
    name_re = re.compile('> \w+')

    # use regex to get names and then to split FASTA
    ma = readall(ma_fasta)
    names = [n.strip('> ') for n in name_re.findall(ma)]
    seqs = [np.array(list(s.replace('\n', ''))) for s in name_re.split(ma)[1:]]

    # check that file parse was correct
    assert len(names) == len(seqs)
    assert all(len(s) == len(seqs[0]) for s in seqs)

    # create dict of name:seq pairs
    fa_dict = dict(izip(names, seqs))

    return fa_dict


def unsplit_alignments(fa_dict, fa_name, nb=100, compress=False):
    """reverse the function 'split alignments', re-save to new FASTA file"""
    # use optional compression
    if compress:
        f = zopen(fa_name, 'w')
    else:
        f = open(fa_name, 'w')
    # write the new file
    for k in fa_dict:
        # reformat and record the sequence header
        f.write('> {}\n'.format(k))
        # split the sequence into nb length lines and write
        seq = fa_dict[k]
        r = range(0, len(seq), nb)
        f.write('\n'.join(''.join(seq[i:i+nb]) for i in r) + '\n')
    f.close()

    return None


def isbase(seq):
    """boolean array checking whether site contains ACGT base"""
    bases = ['A', 'C', 'G', 'T']
    return np.in1d(seq, bases)


def neutmask_slice(cst, pctmin, pctmax, filt=False, all_filters=True):
    """get neutral mask for a range of neutral percentiles"""
    # set file paths
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    fwigz = '{p}/{ncon}/{ncon}.{chrom}.wig.gz'.format(p=pth, **cst.dict)
    if filt:
        fmt = '{p}/{cn}/{cn}.filtered.scores.dist.txt'
    else:
        fmt = '{p}/{cn}/{cn}.scores.dist.txt'
    fdist = fmt.format(p=pth, cn=cst.ncon)
    # get pmin/pmax values based on percentile min/max
    pmin = dtl.dist_percentile(fdist, pctmin)
    pmax = dtl.dist_percentile(fdist, pctmax)

    # get neutral mask using pmin/pmax params
    ncmask = dtl.conservation_mask(fwigz, cst.chrom, pmin=pmin, pmax=pmax)

    # by default, run all filters on the mask
    if all_filters:
        # strict calling mask from 1000 genomes project
        callmask = dtl.sequence_mask(cst.call_mask, cst.chrom)

        # aggregate filter regions into a binary mask array
        filtmask = dtl.get_filter_mask(cst.chrom, cst.files.fnff)

        ncmask &= callmask
        ncmask &= filtmask

    return ncmask


def call_phylofit(ch, multi, pmin, pmax):
    """run phyloFit on neutral sub-alignments of data"""
    # # set alignment label
    # multi = 'primate'

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    # chr1.euarchontoglires.neut_align.0.20
    f_fa = '{}/{}.{}.neut_align.{}.{}.fa'.format(fa_dir, ch, multi, pmin, pmax)
    save_dir = '{}/runs/{}_neutral'.format(phast_dir, multi)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # set phyloFit command parameters
    tree = '{}/{}.tree'.format(cfg_dir, multi)
    pref = '{}/{}.neut_align.{}.{}'.format(save_dir, ch, pmin, pmax)
    smod = 'U2S'
    flag = '-E -p MED -Z -q'

    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=f_fa, op=pref, sm=smod, fl=flag)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)


def call_phylofit_sliding(ch, size=5000, shift=1000):
    """run phyloFit on neutral sub-alignments of data"""
    # record the start time for runtime logging
    start_time = dtime.now()

    # set alignment label
    multi = 'primate'

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    f_fa = '{}/{}.{}.neutral.aligned.4.fa'.format(fa_dir, ch, multi)
    save_dir = '{}/runs/primate_neutral_sliding'.format(phast_dir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # set phyloFit command parameters
    tree = '{}/{}.4.tree'.format(cfg_dir, multi)
    pref = '{}/{}.neutral_sliding'.format(save_dir, ch)
    smod = 'U2S'
    flag = '-E -p MED -Z'
    wind = '-w {},{}'.format(size, shift)

    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl} {wn}'
    phy_cmd = fmt.format(tr=tree, fa=f_fa, op=pref, sm=smod, fl=flag, wn=wind)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))

    return None


def neutral_subrate_percentiles(ch, ncon, pct_min, pct_max):
    """rewrite FASTA and mask files containing only aligned, neutral bases"""
    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_file = pdir + '/fa/{sp}/{c}.{sp}.fa.gz'.format(sp=ncon, c=ch)
    ms_file = pdir + '/fa/{sp}/pct_masks/{c}.{sp}.{pi}.{pj}.nmsk.npz'.\
        format(sp=ncon, c=ch, pi=int(100*pct_min), pj=int(100*pct_max))
    fa_newfile = pdir + '/fa/{sp}/{c}.{sp}.neut_align.{pi}.{pj}.fa'.\
        format(sp=ncon, c=ch, pi=int(100*pct_min), pj=int(100*pct_max))

    # set cst with basic ncons token for loading mask files
    cst = ChromStruct(chrom=ch, ncon=ncon)

    # generate neutral mask for the percentile range
    nmsk = neutmask_slice(cst, pct_min, pct_max)

    # get dict of split alignments as string arrays
    fa_dict = split_alignments(fa_file)

    # get dict of isbase masks for each alignment
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())

    # get the intersection of all isbase masks
    imsk = isbase_dict['hg19']
    for sp in isbase_dict.keys():
        imsk &= isbase_dict[sp]

    # intersect the alignment mask with the neutral mask
    nmsk &= imsk

    # save a copy of the intersect+neutrality mask
    np.savez_compressed(ms_file, neutmask=nmsk)

    # count the initial and final number of sites before and after masking
    init, final = len(nmsk), np.sum(nmsk)
    assert final < init
    msg = 'neutral aligned = {:.4e}'.format(final)
    err_msg(msg)

    # remove the human sequence from set of sequences
    del fa_dict['hg19']

    # keep only neutral sites from the alignment intersection in each seq
    for k in fa_dict:
        fa_dict[k] = fa_dict[k][nmsk]
    # check that the numbers are correct
    assert all(len(v) == final for v in fa_dict.values())

    # re-write the FASTA files
    unsplit_alignments(fa_dict, fa_newfile)


def hco_neutral_phylofit(ch):
    """rewrite FASTA and mask files containing only HCO aligned, neutral"""
    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    # fa_file = pdir + '/fa/primate/{}.primate.fa'.format(ch)
    fa_file = pdir + '/fa/primate/{}.hco.neut.fa'.format(ch)

    # # set cst with basic ncons token for loading mask files
    # cst = ChromStruct(chrom=ch, ncon=ncon)
    #
    # # generate neutral mask for the percentile range
    # nmsk = neutmask_slice(cst, pct_min, pct_max)
    # n_init = nmsk.sum()
    #
    # # get dict of split alignments as string arrays
    # fa_dict = split_alignments(fa_file)
    #
    # # keep only HCO species in fa_dict
    # spc = ['hg19', 'panTro4', 'ponAbe2']
    # fa_dict = dict((s, fa_dict[s]) for s in spc)
    # assert len(fa_dict) == 3
    #
    # # get dict of isbase masks for each alignment
    # isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())
    #
    # # intersect neutmask with isbase masks
    # for sp in isbase_dict.keys():
    #     nmsk &= isbase_dict[sp]
    #
    # # count the initial and final number of sites before and after masking
    # n_final = np.sum(nmsk)
    # assert n_final < n_init
    # msg = 'initial/final {} {}'.format(n_init, n_final)
    # err_msg(msg)
    #
    # # keep only neutral sites from the alignment intersection in each seq
    # for k in fa_dict:
    #     fa_dict[k] = fa_dict[k][nmsk]
    # # check that the numbers are correct
    # assert all(len(v) == n_final for v in fa_dict.values())
    #
    # # re-write the FASTA files
    # unsplit_alignments(fa_dict, fa_newfile)

    # run phylofit on the new FASTA multiple alignment file
    save_dir = '{}/runs/hco_neutral'.format(pdir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # set phyloFit command parameters
    tree = pdir + '/cfgs/HCO.tree'
    pref = '{}/{}.hco'.format(save_dir, ch)
    smod = 'U2S'
    flag = '-E -p MED -Z -q'

    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=fa_file, op=pref, sm=smod, fl=flag)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)


def hco_conserved_phylofit(ch, cons, pct_cons, filt=False):
    """rewrite FASTA and mask files containing only HCO aligned, neutral"""
    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_file = pdir + '/fa/primate/{}.primate.fa'.format(ch)
    fmt = '{}/fa/primate/{}.hco.{}.cons.{}.fa'
    fa_newfile = fmt.format(pdir, ch, cons, pct_cons)
    if filt:
        fa_newfile = fa_newfile.replace('.fa', '.filt.fa')

    # generate conserved mask for the given percentile and species
    cst = ChromStruct(chrom=ch, ncon=cons)
    nmsk = neutmask_slice(cst, pct_cons, 1.0, filt=filt, all_filters=False)
    n_init = nmsk.sum()

    # filter any overlap with neutral if filt flag is used
    if filt:
        cst = ChromStruct(chrom=ch, ncon='euarchontoglires')
        fmsk = neutmask_slice(cst, 0.0, 0.3, all_filters=False)
        nmsk &= ~fmsk
        if cons == 'euarchontoglires':
            assert nmsk.sum() == n_init
        n_init = nmsk.sum()

    # get dict of split alignments as string arrays
    fa_dict = split_alignments(fa_file)

    # keep only HCO species in fa_dict
    spc = ['hg19', 'panTro4', 'ponAbe2']
    fa_dict = dict((s, fa_dict[s]) for s in spc)
    assert len(fa_dict) == 3

    # get dict of isbase masks for each alignment
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())

    # intersect neutmask with isbase masks
    for sp in isbase_dict.keys():
        nmsk &= isbase_dict[sp]

    # count the initial and final number of sites before and after masking
    n_final = np.sum(nmsk)
    assert n_final < n_init
    msg = 'initial/final {} {}'.format(n_init, n_final)
    err_msg(msg)

    # keep only neutral sites from the alignment intersection in each seq
    for k in fa_dict:
        fa_dict[k] = fa_dict[k][nmsk]
    # check that the numbers are correct
    assert all(len(v) == n_final for v in fa_dict.values())

    # re-write the FASTA files
    unsplit_alignments(fa_dict, fa_newfile)

    # run phylofit on the new FASTA multiple alignment file
    save_dir = '{}/runs/hco_neutral'.format(pdir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # set phyloFit command parameters
    tree = pdir + '/cfgs/HCO.tree'
    if filt:
        pref = '{}/{}.hco.{}.cons.{}.filt'.format(save_dir, ch, cons, pct_cons)
    else:
        pref = '{}/{}.hco.{}.cons.{}'.format(save_dir, ch, cons, pct_cons)
    smod = 'U2S'
    flag = '-E -p MED -Z -q'

    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=fa_file, op=pref, sm=smod, fl=flag)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)


def hco_reruns(ch, cons, pct_cons, filt=False):
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fmt = '{}/fa/primate/{}.hco.{}.cons.{}.fa'
    fa_file = fmt.format(pdir, ch, cons, pct_cons)
    if filt:
        fa_file = fa_file.replace('.fa', '.filt.fa')
    # run phylofit on the new FASTA multiple alignment file
    save_dir = '{}/runs/hco_neutral'.format(pdir)

    # set phyloFit command parameters
    tree = pdir + '/cfgs/HCO.tree'
    if filt:
        pref = '{}/{}.hco.{}.cons.{}.filt'.format(save_dir, ch, cons, pct_cons)
    else:
        pref = '{}/{}.hco.{}.cons.{}'.format(save_dir, ch, cons, pct_cons)
    smod = 'U2S'
    flag = '-E -p MED -Z -q'

    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=fa_file, op=pref, sm=smod, fl=flag)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)


def conserved_annotation_overlaps(ch, pct_cons, filt=False):
    """rewrite FASTA and mask files containing only HCO aligned, neutral"""
    # generate conserved mask for the given percentile and species
    msks = {}
    cons = ['ape', 'euarchontoglires', 'fish']
    for cn in cons:
        cst = ChromStruct(chrom=ch, ncon=cn)
        nmsk = neutmask_slice(cst, pct_cons, 1.0, filt, False)
        n_init = nmsk.sum()
        # filter any overlap with neutral if filt flag is used
        if filt:
            cst = ChromStruct(chrom=ch, ncon='euarchontoglires')
            fmsk = neutmask_slice(cst, 0.0, 0.3, all_filters=False)
            nmsk &= ~fmsk
            if cons == 'euarchontoglires':
                assert nmsk.sum() == n_init
            # n_init = nmsk.sum()
        msks[cn] = nmsk
    for cn in cons:
        err_msg('{} {}'.format(cn, msks[cn].sum()))
    err_msg('A+E {}'.format(np.sum(msks['ape'] & msks['euarchontoglires'])))
    err_msg('A+F {}'.format(np.sum(msks['ape'] & msks['fish'])))
    err_msg('E+F {}'.format(np.sum(msks['euarchontoglires'] & msks['fish'])))
    aef = np.sum(msks['euarchontoglires'] & msks['fish'] & msks['ape'])
    err_msg('A+E+F {}'.format(aef))


def main_4():
    """rewrite FASTA and mask files containing only aligned, neutral bases"""
    # get chrom from command line
    if len(argv) != 2:
        print 'usage: align_stats <ch>'
        exit(1)
    ch = argv[1]

    # set FASTA file dir
    fa_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/primate'
    # set cst with basic ncons token for loading mask files
    cst = ChromStruct(chrom=ch, tkn='basic.ncons')
    # load neutmask
    nmsk = np.load(cst.neut_masks)['neutmask']

    # set the list of species used in the new alignment
    # spc = ['gorGor3', 'ponAbe2', 'panTro4', 'papHam1', 'chlSab1', 'macFas5',
    #        'nomLeu3', 'rheMac3']
    spc = ['panTro4', 'chlSab1', 'macFas5']
    # hg19_panTro4_chlSab1_macFas5
    # get dict of split alignments as string arrays
    fa_file = '{}/{}.primate.fa'.format(fa_dir, ch)
    fa_dict = split_alignments(fa_file)

    # get data from each species and a selection of species intersections
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())

    # get the intersection of all isbase masks
    imsk = isbase_dict['hg19']
    for sp in spc:
        imsk &= isbase_dict[sp]

    # intersect the alignment mask with the neutral mask
    nmsk &= imsk

    # save a copy of the intersect+neutrality mask
    new_tkn = 'basic.ncons.aligned.{}'.format(len(spc)+1)
    fmask = cst.neut_masks.replace('basic.ncons', new_tkn)
    np.savez(fmask, neutmask=nmsk)

    # count the initial and final number of sites before and after masking
    init, final = len(nmsk), np.sum(nmsk)
    assert final < init

    # keep only neutral sites from the alignment intersection in each seq
    for k in fa_dict:
        fa_dict[k] = fa_dict[k][nmsk]
    # check that the numbers are correct
    assert all(len(v) == final for v in fa_dict.values())

    # re-write the FASTA files
    fa_tkn = '.neutral.aligned.{}.fa'.format(len(spc)+1)
    fa_newfile = fa_file.replace('.fa', fa_tkn)
    unsplit_alignments(fa_dict, fa_newfile)


def main_3():
    """call phyloFit using sliding winows"""
    # get chrom from command line
    if len(argv) != 2:
        print 'usage: align_stats <ch>'
        exit(1)
    ch = argv[1]
    size = 100000
    shift = 20000
    call_phylofit_sliding(ch, size, shift)


def main_1():
    """rewrite FASTA alignment files containing only aligned, neutral bases"""
    # get chrom from command line
    if len(argv) != 2:
        print 'usage: align_stats <ch>'
        exit(1)
    ch = argv[1]

    # set FASTA file dir
    fa_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/primate'
    # set cst with basic ncons token for loading mask files
    cst = ChromStruct(chrom=ch, tkn='basic.ncons')
    # load neutmask
    nmsk = np.load(cst.neut_masks)['neutmask']

    # get dict of split alignments as string arrays
    fa_file = '{}/{}.primate.fa'.format(fa_dir, ch)
    fa_dict = split_alignments(fa_file)

    # get data from each species and a selection of species intersections
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())

    # get the intersection of all isbase masks
    imsk = isbase_dict['hg19']
    spc = ['gorGor3', 'ponAbe2', 'panTro4', 'papHam1', 'chlSab1', 'macFas5',
           'nomLeu3', 'rheMac3']
    for sp in spc:
        imsk &= isbase_dict[sp]

    # intersect the alignment mask with the neutral mask
    nmsk &= imsk

    # count the initial and final number of sites before and after masking
    init, final = len(nmsk), np.sum(nmsk)
    assert final < init

    # keep only neutral sites from the alignment intersection in each seq
    for k in fa_dict:
        fa_dict[k] = fa_dict[k][nmsk]
        # fa_dict[k] = fa_dict[k][nmsk[imsk]]
    # check that the numbers are correct
    # assert all(len(v) == final for v in fa_dict.values())
    assert all(len(v) == final for v in fa_dict.values())

    # re-write the FASTA files
    fa_newfile = fa_file.replace('.fa', '.neutral.fa')
    unsplit_alignments(fa_dict, fa_newfile)

    # record the initial and final number of sites in the new FASTA files
    flog = '{}/{}.rewrite.log'.format(fa_dir, ch)
    with open(flog, 'w') as f:
         f.write('BASE COUNT: INITIAL={} FINAL={}\n'.format(init, final))


def main_0():
    """print percentage alignment for different combinations of species"""
    # get chrom from command line
    if len(argv) != 2:
        print 'usage: align_stats <ch>'
        exit(1)
    ch = argv[1]

    # set FASTA file dir
    fa_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/primate'
    # set cst with basic ncons token for loading mask files
    cst = ChromStruct(chrom=ch, tkn='basic.ncons')
    # load neutmask
    nmsk = np.load(cst.neut_masks)['neutmask']

    # get dict of split alignments as string arrays
    fa_file = '{}/{}.primate.fa'.format(fa_dir, ch)
    fa_dict = split_alignments(fa_file)

    # get data from each species and a selection of species intersections
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())
    # get human macaque intersection
    hm_itx = isbase_dict['hg19'] & isbase_dict['rheMac3']

    # # record initial all human as well as HM values
    # labels = ['n', 'hm']
    # values = [np.sum(nmsk[isbase_dict['hg19']]), hm]

    # record counts for each alignment configuration to a file
    fout = '/ifs/data/c2b2/gs_lab/dam2214/run/{}.align_stats.log'.format(ch)
    with open(fout, 'w') as f:
        # record initial all human as well as HM values
        f.write('n {}\n'.format(np.sum(nmsk[isbase_dict['hg19']])))
        f.write('hm {}\n'.format(np.sum(nmsk[hm_itx])))

        # create check the remaining bases for every alignment configuration
        spc = ['gorGor3', 'ponAbe2', 'panTro4', 'papHam1', 'chlSab1', 'macFas5',
               'nomLeu3', 'rheMac3']
        for n in range(1, len(spc)+1):
            for cmb in combinations(spc, n):
                l = 'hg19'
                h = np.copy(isbase_dict['hg19'])
                for sp in cmb:
                    l += '_{}'.format(sp)
                    h &= isbase_dict[sp]
                f.write('{} {}\n'.format(l, np.sum(nmsk[h])))
                # labels.append(l)
                # values.append(np.sum(nmsk[h]))


if __name__ == '__main__':
    # if len(argv) != 4:
    #     print 'usage: align_stats <chrom> <pctmin> <pctmax>'
    #     exit(1)
    # pctmin, pctmax = map(float, argv[2:])
    # call_phylofit(ch, ncon, pctmin, pctmax)

    # if len(argv) != 5:
    #     print 'usage: align_stats <chrom> <cons> <pct_cons> <filter>'
    #     exit(1)
    # ch = argv[1]
    # cons = argv[2]
    # pct_cons = float(argv[3])
    # filt = eval(argv[4])
    # hco_reruns(ch, cons, pct_cons, filt)

    if len(argv) != 4:
        print 'usage: align_stats <chrom> <pct_cons> <filter>'
        exit(1)
    ch = argv[1]
    pct_cons = float(argv[2])
    filt = eval(argv[3])
    conserved_annotation_overlaps(ch, pct_cons, filt)

    # if len(argv) != 2:
    #     print 'usage: align_stats <chrom>'
    #     exit(1)
    # ch = argv[1]
    # hco_neutral_phylofit(ch)