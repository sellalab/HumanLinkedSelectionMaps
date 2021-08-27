__author__ = 'davidmurphy'


import os
import re
import subprocess
import numpy as np
from sys import stderr, stdout, argv
from gzip import open as zopen
from collections import Counter
from datetime import datetime as dtime
from itertools import izip, combinations
from classes.runstruct import ChromStruct
import data_processing.data_tools as dtl
from data_processing.parse_cadd import CADDScore


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


def isbase(seq):
    """boolean array checking whether site contains ACGT base"""
    bases = ['A', 'C', 'G', 'T']
    return np.in1d(seq, bases)


def npad_alignment(fa, msk):
    """replace gaps with 2 Ns in alignment array"""
    # loop through segments and add two Ns to the start of each one
    segments = []
    for (start, end) in dtl.binary_mask_segments(msk):
        # create the padded segment
        pad_seg = np.concatenate((np.array(['N', 'N']), fa[start:end]))
        segments.append(pad_seg)

    # convert segments into a single array
    segments = np.concatenate(segments)

    return segments


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


def join_alignments(fa_dict, fa_name, nb=100, compress=False):
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
        r = xrange(0, len(seq), nb)
        f.write('\n'.join(''.join(seq[i:i+nb]) for i in r) + '\n')
    f.close()

    return None


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


def parse_exptotsub(xts_file, n_cells=16):
    """parse the results of expected total substitution matrices for MA"""
    # read lines into a single list, cut off the header
    with open(xts_file, 'r') as f:
        lines = f.readlines()[8:]

    # save matrix to dict hashed by node label
    mat_dict = {}

    # create temp variables for each matrix
    cur_mat = []
    cur_lbl = ''

    # extract each matrix and its label from the lines
    for line in lines:
        # split line on whitespace
        line = line.split()

        # skip blank lines
        if len(line) == 0:
            continue

        # get node number and name
        if line[0] == 'Branch':
            cur_lbl = line[-1].split("'")[1]

        # record the matrix
        if len(line) == n_cells + 1:
            cur_mat.append(map(float, line[1:]))

        # when current 16x16 matrix is complete, store in dict
        if len(cur_mat) == n_cells:
            mat_dict[cur_lbl] = np.array(cur_mat)
            cur_mat = []

    return mat_dict


def hcgo_neutral_subalignment(ch):
    """build a subalignment of neutral segments (padded with NN) for analysis"""
    # fix some variables
    ncon = 'euarchontoglires'
    pct_min, pct_max = 0, 0.35

    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_file = pdir + '/fa/primate/{}.primate.fa'.format(ch)
    fa_subalign = pdir + '/fa/primate/{}.hcgo.neut.fa'.format(ch)

    # set cst with basic ncons token for loading mask files
    cst = ChromStruct(chrom=ch, ncon=ncon)

    # generate neutral mask for the percentile range
    nmsk = neutmask_slice(cst, pct_min, pct_max, filt=True)
    n_init = nmsk.sum()

    # get dict of split alignments as string arrays
    fa_dict = split_alignments(fa_file)

    # keep only HCGO species in fa_dict
    spc = ['hg19', 'panTro4', 'gorGor3', 'ponAbe2']
    fa_dict = dict((s, fa_dict[s]) for s in spc)
    assert len(fa_dict) == 4

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
    for (k, fa) in fa_dict.items():
        # add N pad to runs of aligned bases for context analysis
        fa_dict[k] = npad_alignment(fa, nmsk)
        # fa_dict[k] = fa_dict[k][nmsk]
    # check that the numbers are correct
    # assert all(len(v) == n_final for v in fa_dict.values())

    # re-write the FASTA files
    join_alignments(fa_dict, fa_subalign)


def hcgo_conserved_subalign(ch, cons, pct_cons, filt=False):
    """rewrite FASTA and mask files containing only HCO aligned, neutral"""
    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_file = pdir + '/fa/primate/{}.primate.fa'.format(ch)
    fmt = '{}/fa/primate/{}.hcgo.{}.cons.{}.fa'
    fa_subalign = fmt.format(pdir, ch, cons, pct_cons)
    if filt:
        fa_subalign = fa_subalign.replace('.fa', '.filt.fa')

    # CADD score rule:
    if cons == 'cadd':
        # use the CADDScore class to get the mask
        cadd_root = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/cadd'
        cadd_dir2 = cadd_root + '/cadd_v1.6_without_bstat'
        score = CADDScore(ch, cadd_dir2)
        cmsk = score.get_mask(float(pct_cons))
    else:
        # USE NEW FORMAT CONSERVATION SCORES
        pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
        # farra = pth + '/{cn}/{cn}.{c}.scores.npz'.format(cn=cons, c=ch)
        # fdist = pth + '/{cn}/{cn}.{c}.clean.score.dist.txt'.format(cn=cons,
        #                                                            c=ch)
        # USE *NEWEST* FORMAT WITH GMAP MASK
        farra = pth + '/{cn}/{cn}.{c}.gmask.scores.npz'.format(cn=cons, c=ch)
        fdist = pth + '/{cn}/{cn}.{c}.gmask.score.dist.txt'.format(cn=cons,
                                                                   c=ch)

        # get the minimum conservation score for the given percentile
        pmin = dtl.dist_percentile(fdist, float(pct_cons) * 0.01)
        # get conservation mask at the given percentage conserved
        cmsk = (np.load(farra)['p'] >= pmin)

    msg = 'conserved_sites {}'.format(np.sum(cmsk))
    err_msg(msg)
    # # convert pct_cons to fraction of 1.0
    # pcons = float(pct_cons) / 100.0
    # nmsk = neutmask_slice(cst, pcons, 1.0, filt=filt, all_filters=False)
    # n_init = nmsk.sum()
    #
    # # filter any overlap with neutral if filt flag is used
    # if filt:
    #     cst = ChromStruct(chrom=ch, ncon='euarchontoglires')
    #     fmsk = neutmask_slice(cst, 0.0, 0.35, all_filters=False)
    #     nmsk &= ~fmsk
    #     if cons == 'euarchontoglires':
    #         assert nmsk.sum() == n_init
    #     n_init = nmsk.sum()

    # get dict of split alignments as string arrays
    fa_dict = split_alignments(fa_file)

    # keep only HCGO species in fa_dict
    spc = ['hg19', 'panTro4', 'gorGor3', 'ponAbe2']
    fa_dict = dict((s, fa_dict[s]) for s in spc)
    assert len(fa_dict) == 4

    # get dict of isbase masks for each alignment
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())

    # intersect neutmask with isbase masks
    for sp in isbase_dict.keys():
        cmsk &= isbase_dict[sp]

    # # count the initial and final number of sites before and after masking
    # n_final = np.sum(cmsk)
    # assert n_final < n_init
    # msg = 'initial/final {} {}'.format(n_init, n_final)

    msg = 'aligned_conserved_sites {}'.format(np.sum(cmsk))
    err_msg(msg)

    # keep only neutral sites from the alignment intersection in each seq
    for (k, fa) in fa_dict.items():
        # add N pad to runs of aligned bases for context analysis
        fa_dict[k] = npad_alignment(fa, cmsk)

    # re-write the FASTA files
    join_alignments(fa_dict, fa_subalign)


def hcgo_neutral_phylofit(ch, mod):
    """rewrite FASTA and mask files containing only HCO aligned, neutral"""
    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_subalign = pdir + '/fa/primate/{}.hcgo.neut.fa'.format(ch)

    # create sub-alignment if it does not yet exist
    if not os.path.isfile(fa_subalign):
        hcgo_neutral_subalignment(ch)

    # run phylofit on the new FASTA multiple alignment file
    save_dir = '{}/runs/hcgo_neutral'.format(pdir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # set phyloFit command parameters
    tree = pdir + '/cfgs/HCGO.tree'
    pref = '{}/{}.hcgo.neut.{}'.format(save_dir, ch, mod)
    smod = mod
    flag = '-E -p MED -Z -q'

    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=fa_subalign, op=pref, sm=smod, fl=flag)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)


def hcgo_conserved_phylofit(ch, cons, pct_cons, mod, filt=False):
    """rewrite FASTA and mask files containing only HCO aligned, neutral"""
    # set file paths
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    ffmt = '{}/fa/primate/{}.hcgo.{}.cons.{}.fa'
    fa_subalign = ffmt.format(pdir, ch, cons, pct_cons)
    if filt:
        fa_subalign = fa_subalign.replace('.fa', '.filt.fa')

    # # create sub-alignment if it does not yet exist
    # if not os.path.isfile(fa_subalign):
    # create sub-alignment, overwrite existing
    msg = 'creating {} {} cons{} sub-alignment'.format(ch, cons, pct_cons)
    err_msg(msg)
    hcgo_conserved_subalign(ch, cons, pct_cons, filt)

    # run phylofit on the new FASTA multiple alignment file
    save_dir = '{}/runs/hcgo_cons'.format(pdir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # set phyloFit command parameters
    tree = pdir + '/cfgs/HCGO.tree'
    ofmt = '{}/{}.hcgo.{}.cons.{}.{}'
    pref = ofmt.format(save_dir, ch, cons, pct_cons, mod)
    if filt:
        pref += '.filt'
    smod = mod
    flag = '-E -p MED -Z -q'

    # run phyloFit on the sub-alignment
    cfmt = 'phyloFit -t {tr} -i FASTA {fa} -o {op} -s {sm} {fl}'
    phy_cmd = cfmt.format(tr=tree, fa=fa_subalign, op=pref, sm=smod, fl=flag)

    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)


def main():
    if len(argv) == 3:
        ch, mod = argv[1:]
        hcgo_neutral_phylofit(ch, mod)
    elif len(argv) == 5:
        ch, cons, pct_cons, mod = argv[1:]
        hcgo_conserved_phylofit(ch, cons, pct_cons, mod)
    else:
        print 'usage_1: neutcons_subrates <chrom> <model>'
        print 'usage_2: neutcons_subrates <chrom> <cons> <pcons> <model>'
        exit(1)


if __name__ == '__main__':
    main()
