__author__ = 'davidmurphy'


import os
import re
import time
import subprocess
import numpy as np
from gzip import open as zopen
from collections import defaultdict
from sys import stderr, stdout, argv
from datetime import datetime as dtime
from classes.phylotree import parse_exptotsub
from classes.runstruct import ChromStruct, chromosome_length, root_dir, \
    human_autosomes, cst_from_fldr


def non_bgc_subs(mat):
    """return a new matric of non-BGC substitution rates"""
    nbgc = 0
    nbgc += mat[0, (3, 12)].sum()  # AA > AT, TA
    nbgc += mat[1, (2, 13)].sum()  # AC > AG, TC
    nbgc += mat[2, (1, 14)].sum()  # AG > AC, TG
    nbgc += mat[3, (0, 15)].sum()  # AT > AA, TT
    nbgc += mat[4, (7, 8)].sum()  # CA > CT, GA
    nbgc += mat[5, (6, 9)].sum()  # CC > CG, GC
    nbgc += mat[6, (5, 10)].sum()  # CG > CC, GG
    nbgc += mat[7, (4, 11)].sum()  # CT > CA, GT
    nbgc += mat[8, (4, 11)].sum()  # GA > CA, GT
    nbgc += mat[9, (5, 10)].sum()  # GC > CC, GG
    nbgc += mat[10, (6, 9)].sum()  # GG > CG, GC
    nbgc += mat[11, (7, 8)].sum()  # GT > CT, GA
    nbgc += mat[12, (0, 15)].sum()  # TA > AA, TT
    nbgc += mat[13, (1, 14)].sum()  # TC > AC, TG
    nbgc += mat[14, (2, 13)].sum()  # TG > AG, TC
    nbgc += mat[14, (3, 12)].sum()  # TT > AT, TA

    return nbgc


def err_msg(errstring):
    stderr.write(errstring + '\n')
    stdout.flush()


def file_aware_reader(f_name):
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


def fix_chr20():
    """rebuild chr20 alignment"""
    spc_8 = ['hg19', 'gorGor3', 'ponAbe2', 'panTro4', 'papHam1', 'chlSab1',
             'macFas5', 'rheMac3']
    f = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/euarchontoglires/chr20.euarchontoglires.fa.gz'
    fa_dict = split_alignments(f)

    # set FASTA file dir
    fa_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/primate'

    # get dict of split alignments as string arrays
    fa_save = '{}/{}.primate.fa'.format(fa_dir, 'chr20')
    new_dict = dict((k, fa_dict[k]) for k in spc_8)
    join_alignments(new_dict, fa_save)

    return None


def split_alignments(ma_fasta):
    """split a multiple alignment FASTA file and return dict of name:seq"""
    # set FASTA label regex
    name_re = re.compile('> \w+')

    # use regex to get names and then to split FASTA
    ma = file_aware_reader(ma_fasta)
    names = [n.strip('> ') for n in name_re.findall(ma)]
    seqs = [np.array(list(s.replace('\n', ''))) for s in name_re.split(ma)[1:]]

    # check that file parse was correct
    assert len(names) == len(seqs)
    assert all(len(s) == len(seqs[0]) for s in seqs)

    # create dict of name:seq pairs
    fa_dict = dict(zip(names, seqs))

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
        r = range(0, len(seq), nb)
        f.write('\n'.join(''.join(seq[i:i+nb]) for i in r) + '\n')
    f.close()

    return None


def process_matrix(fname, return_string=True):
    """get sub counts from exptotsub matrix file and return formatted line"""
    # get chrom and base range by splitting filename
    ch, rng = fname.split('/')[-1].split('.')[:2]
    # get start, end positions from the range
    start, end = map(int, rng.split('_'))

    # process the expected substitution matrices from file
    dmat = parse_exptotsub(fname)
    # get the matrices for ALL branches in the tree
    mats = dmat.values()

    # calculate the total number of bases across branches
    bases = sum(m.sum() for m in mats)
    # calculate total substitutions across branches
    subs = bases - sum(m.diagonal().sum() for m in mats)
    # calcuate CpG substitutions across branches
    cpg = sum(m[6, (4, 14)].sum() for m in mats)
    # calculate non-BGC substitutions across branches
    nbgc = sum(non_bgc_subs(m) for m in mats)

    # adjust bases by number of trees
    try:
        bases /= len(mats)
    except ZeroDivisionError:
        err_msg('MATRIX FILE {} TRIGGERED ZERO DIVISION ERROR'.format(fname))

    # count total bases, subs and CpG adjusted subs
    b_cnt = int(bases)
    s_cnt = int(0.5 * subs)
    c_cnt = int(0.5 * (subs - cpg))
    g_cnt = int(0.5 * nbgc)

    # format results from current range to append to file
    fmt_1 = '{st} {en} {bs} {sb}\n'
    fmt_2 = '{st} {en} {bs} {sb} {cp} {bg}\n'

    line_1 = fmt_1.format(st=start, en=end, bs=b_cnt, sb=s_cnt)
    line_2 = fmt_2.format(st=start, en=end, bs=b_cnt, sb=s_cnt, cp=c_cnt,
                          bg=g_cnt)

    # return data string
    if return_string:
        return line_2
    # return tuple of data
    else:
        return start, end, b_cnt, s_cnt


def create_neutral_alignment(ch, keep_species):
    """rewrite FASTA and mask files containing only aligned, neutral bases"""
    # set FASTA file dir
    fa_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/primate'
    # use a neutmask with the euarchontoglires 30% conserved sites
    cst = ChromStruct(chrom=ch)

    # load neutmask (get initial sum of passing sites)
    nmsk = np.load(cst.neut_masks)['neutmask']
    init = np.sum(nmsk)

    # get dict of split alignments as string arrays
    fa_file = '{}/{}.primate.fa'.format(fa_dir, ch)
    fa_dict = split_alignments(fa_file)

    # keep only the species of interest in the dictionary
    fa_dict = dict((s, fa_dict[s]) for s in keep_species)

    # get data from each species and a selection of species intersections
    isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())

    # get the intersection of all isbase masks with nmsk
    for sp in keep_species:
        nmsk &= isbase_dict[sp]
    # (get final sum of passing sites)
    final = np.sum(nmsk)
    # count the initial and final number of sites before and after masking
    assert final < init
    msg = 'initial {} final {}\n'.format(init, final)
    stderr.write(msg)
    stdout.flush()

    # save a copy of the intersect+neutrality mask
    new_tkn = '.aligned.{}.nmsk'.format(len(keep_species))
    f_mask = cst.neut_masks.replace('.nmsk', new_tkn)
    np.savez_compressed(f_mask, neutmask=nmsk)

    # keep only neutral sites from the alignment intersection in each seq
    for k in fa_dict:
        fa_dict[k] = fa_dict[k][nmsk]
    # check that the numbers are correct
    assert all(len(v) == final for v in fa_dict.values())

    # re-write the FASTA files
    fa_tkn = '.neutral.aligned.{}.fa'.format(len(keep_species))
    fa_newfile = fa_file.replace('.fa', fa_tkn)
    join_alignments(fa_dict, fa_newfile)


def create_bbin_alignment(cst, bbin):
    """create alignment of all sites for a given B bin"""
    # record start time for elapsed time calculations
    start_time = dtime.now()
    # kept species for the pseudo-alignment
    spc_8 = ['hg19', 'gorGor3', 'ponAbe2', 'panTro4', 'papHam1', 'chlSab1',
             'macFas5', 'rheMac3']
    # set FASTA dir
    fa_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/fa/primate'
    # set save file
    fa_save = '{}/primate.{}.bbin{}.fa'.format(fa_dir, cst.tkn, bbin)

    # load data for the B bin
    bin_dir = root_dir + '/compress/sorted_bins'
    f_load = bin_dir + '/sorted_bin{}.{}.npz'.format(bbin, cst.tkn)
    bdata = np.load(f_load)['sg'].astype(int)

    # build the pseudo-alignment for the B bin one chrom at a time
    bin_site_count = 0
    masked_bin_site_count = 0
    bin_dict = defaultdict(list)
    for ch in human_autosomes:
        # update chromosome on the cst
        cst.chrom = ch
        # load the neutral mask to get site positions
        nmsk = np.load(cst.neut_masks)['neutmask']

        # get mask of all segments for this chrom
        ch_int = int(ch[3:])
        cmsk = (bdata[:,0] == ch_int)
        # continue if there is no data from the current chrom
        if np.sum(cmsk) == 0:
            continue

        # string together all of the neutral positions for the current chrom
        pos = []
        for (start, end, num) in bdata[cmsk][:, 1:]:
            npos = np.where(nmsk[start:end] > 0)[0] + start
            assert npos.size == num
            pos.append(npos)
        # concatenate all of the neutral positions for the current chrom
        pos = np.concatenate(pos)
        # track the sites in this bin
        c_sites = len(pos)
        bin_site_count += c_sites

        # get dict of split alignments as string arrays
        fa_file = '{}/{}.primate.fa'.format(fa_dir, ch)
        fa_dict = split_alignments(fa_file)
        # keep only the species of interest in the dictionary
        fa_dict = dict((s, fa_dict[s]) for s in spc_8)
        # keep only the neutral sites for the current chrom in the dictionary
        for k in fa_dict:
            fa_dict[k] = fa_dict[k][pos]
        # create a mask that weeds out non-bases from the new alignment
        isbase_dict = dict((k, isbase(fa_dict[k])) for k in fa_dict.keys())
        bmsk = isbase_dict[spc_8[0]]
        assert len(bmsk) == pos.size
        for k in spc_8[1:]:
            bmsk &= isbase_dict[k]
        # track the number of sites after masking
        mc_sites = np.sum(bmsk)
        masked_bin_site_count += mc_sites

        # add masked arrays for current chrom to the master bin dictionary
        for k in spc_8:
            bin_dict[k].append(fa_dict[k][bmsk])

        # print status message for the chromosome just completed
        msg = '{} neut {} masked {}'.format(ch, c_sites, mc_sites)
        err_msg(msg)

    # concatenate the positions from each of the chromosomes for each seq
    for k in spc_8:
        bin_dict[k] = np.concatenate(bin_dict[k])
    # join sequences in a new pseudo-alignment
    join_alignments(bin_dict, fa_save)

    # run phyloFit on the pseudo-alignment
    process_bbin_subalign(cst, bbin)

    # print run stats
    msg = 'neut {} masked {}'.format(bin_site_count, masked_bin_site_count)
    err_msg(msg)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))

    return None


def check_bbin():
    """create alignment of all sites for a given B bin"""
    # load data for the B bin
    bin_dir = root_dir + '/compress/sorted_bins'
    # f_save = bin_dir + '/{}.pos.npz'
    f_save = bin_dir + '/{}.seg.npz'

    # for bbin in range(100):
    #     bbin_pos = []
    #     f_load = bin_dir + '/sorted_bin{}.npz'.format(bbin)
    #     bdata = np.load(f_load)['sg']
    #     err_msg('bbin {} seglen {}'.format(bbin, len(bdata)))

    # reconstruct all of th bbin segments for each chrom
    for ch in human_autosomes:
        seg_len = 0
        pos = []
        seg = []
        # load the neutral mask to get site positions
        cst = ChromStruct(chrom=ch)
        nmsk = np.load(cst.neut_masks)['neutmask']
        # get mask of all segments for this chrom
        ch_int = int(ch[3:])
        for bbin in range(100):
            bbin_pos = []
            f_load = bin_dir + '/sorted_bin{}.npz'.format(bbin)
            bdata = np.load(f_load)['sg']
            cmsk = (bdata[:,0] == ch_int)
            # continue if there is no data from the current chrom
            if np.sum(cmsk) == 0:
                continue
            # string together all of the neutral positions for the current chrom
            for (start, end) in bdata[cmsk][:, 1:]:
                seg_len += (end-start)
                npos = np.where(nmsk[start:end] > 0)[0] + start
                bbin_pos.append(npos)
            # concatenate all of the neutral positions for the current chrom
            bbin_pos = np.concatenate(bbin_pos)
            # add bbin pos to total pos list for chrom
            pos.append(bbin_pos)
            seg.append(bdata[cmsk][:, 1:])

        seg = np.concatenate(seg)
        # sort segs by first position
        sidx = np.argsort(seg[:,0])

        # sort pos and save
        # pos = np.concatenate(pos).astype('u4')
        # pos.sort()

        # np.savez_compressed(f_save.format(ch), seg=seg[sidx])
        err_msg('{}: recorded={} neut={} seglen={}'.format(ch, len(pos), np.sum(nmsk), seg_len))

    return None


def process_bbin_subalign(cst, bbin):
    """process the B bin pseudo-alignments with phyloFit"""
    # # record the start time for runtime logging
    # start_time = dtime.now()

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = phast_dir + '/fa/primate'
    save_dir = phast_dir + '/runs/bbins'
    # create save dir if needed
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = phast_dir + '/cfgs'

    # set pseudo-alignment file path
    fa_load = '{}/primate.{}.bbin{}.fa'.format(fa_dir, cst.tkn, bbin)
    # set phyloFit command parameters
    tree = cfg_dir + '/primate.8.tree'
    # use physical positions for save label preference
    pref = '{}/{}.bbin{}'.format(save_dir, cst.tkn, bbin)
    smod = 'UNREST'
    flag = '-E -p MED -Z -q'
    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {lb} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=fa_load, lb=pref, sm=smod, fl=flag)
    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)

    # # calculate the runtime and write to log
    # run_time = dtime.now() - start_time
    # err_msg('RUN TIME: {}'.format(run_time))

    return None


def bbin_subcounts():
    """get AT and GC subcounts for each bin"""
    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    save_dir = phast_dir + '/runs/bbins'
    res = []
    for bbin in range(100):
        # get matrix of sub counts
        mat = np.zeros(shape=(4, 4))
        f_name = '{}/bbin{}.exptotsub'.format(save_dir, bbin)
        dmat = parse_exptotsub(f_name, n_cells=4)
        # add all branches into one matrix
        for m in dmat.values():
            mat += m

        # get AT substitution rate
        at = sum(mat[0, 1:]) + sum(mat[3, 0:3])
        at_tot = np.sum(mat[(0, 3), :]) / len(dmat.values())
        at_rate = at / at_tot

        # get ONLY A>T mutations as a special subset
        a2t_rate = mat[0, 3] * len(dmat.values()) / np.sum(mat[0, :])

        # get GC substitution rate
        gc = sum(mat[1, (0, 2, 3)]) + sum(mat[2, (0, 1, 3)])
        gc_tot = np.sum(mat[1:3, :]) / len(dmat.values())
        gc_rate = gc / gc_tot

        # get substitution rate for all mutations
        all_rate = (at + gc) / (at_tot + gc_tot)

        # get HM substitution rate
        hm_br = [k for k in dmat.keys() if ('hg19' in k) or ('rheMac3' in k)]
        hm_mat = np.zeros(shape=(4, 4))
        for k in hm_br:
            hm_mat += dmat[k]
        hm = np.sum(hm_mat) - np.sum(np.diagonal(hm_mat))
        hm_tot =  np.sum(hm_mat) / len(hm_br)
        hm_rate = hm / hm_tot

        # create row of all rates
        row = (bbin, at_rate, a2t_rate, gc_rate, all_rate, hm_rate)
        res.append(row)

    # convert results to an array and save
    res = np.array(res)
    f_save = save_dir + '/bbin_rates.txt'
    np.savetxt(f_save, res, fmt='%d %f %f %f %f %f')

    return None


def process_neutral_segment(ch, start, end, win):
    # record the start time for runtime logging
    start_time = dtime.now()

    # fix "multi" to "primate"
    multi = 'primate'
    nspec = 8

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = phast_dir + '/fa/primate'
    save_dir = '{}/runs/aligned{}_win{}'.format(phast_dir, nspec, win)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # format output file name for substitution counts file
    fout = '{}/{}.{}_{}.subcnt.txt'.format(save_dir, ch, start, end)

    # set path to input FASTA file
    fa_in = '{}/{}.primate.neutral.aligned.{}.fa'.format(fa_dir, ch, nspec)

    # use a neutmask with the euarchontoglires 30% conserved sites
    cst = ChromStruct(chrom=ch)
    # nmsk = np.load(cst.neut_masks)['neutmask']
    new_tkn = '.aligned.{}.nmsk'.format(nspec)
    f_mask = cst.neut_masks.replace('.nmsk', new_tkn)
    nmsk = np.load(f_mask)['neutmask']

    # get positions of neutral sites from the neutral mask
    npos = np.where(nmsk == 1)[0]

    # for the start end range, create 20kb sub-alignments for phyloFit
    for i in range(start, end, win):
        # set upper bound for region; don't exceed total length of FASTA
        j = min(len(npos) - 1, i + win)
        # get PHYSICAL positions of i & j for writing subcount file
        ipos, jpos = npos[i], npos[j]

        # set output sub-alignment FASTA file name
        fa_sub = fa_in.replace('.fa', '.{}_{}.fa'.format(i + 1, j))
        # create W-kb sub-alignment using msa_view program
        fmt = 'msa_view {fin} --start {st} --end {en} --refidx 1 > {fout}'
        msa_cmd = fmt.format(fin=fa_in, st=i + 1, en=j, fout=fa_sub)
        # call msa_view through the shell
        subprocess.call(msa_cmd, shell=True)

        # set phyloFit command parameters
        tree = '{}/primate.{}.tree'.format(cfg_dir, nspec)
        # use physical positions for save label preference
        pref = '{}/{}.{}_{}'.format(save_dir, ch, ipos, jpos)
        smod = 'U2S'
        flag = '-E -p MED -Z -q'
        # run phyloFit on the sub-alignment
        fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {lb} -s {sm} {fl}'
        phy_cmd = fmt.format(tr=tree, fa=fa_sub, lb=pref, sm=smod, fl=flag)
        # call phyloFit through the shell
        subprocess.call(phy_cmd, shell=True)

        # set phyloFit output file names
        mat_file = '{}.exptotsub'.format(pref)
        mod_file = '{}.mod'.format(pref)

        # use try/except to continue through IOErrors
        try:
            # process phyloFit file to get sub counts and append to out file
            scnt_line = process_matrix(mat_file)
            with open(fout, 'a') as f:
                f.write(scnt_line)
            # pause 1 second before file deletion because of unsync errors
            time.sleep(2)
            # delete unused mod file
            os.remove(mod_file)
            # delete matrix file after processing
            os.remove(mat_file)
            # delete the sub-alignment
            os.remove(fa_sub)

        except IOError:
            # err_msg('WARNING: {} missing.'.format(mat_file))
            err_msg('WARNING: {} {} {} MISSING'.format(ch, i, j))
        except OSError:
            err_msg('WARNING: {} {} {} MISSING'.format(ch, i, j))

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))


def combine_subcount_files(ch, win):
    # set directories used in run
    nspec = 8
    pdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    pth = '{}/runs/aligned{}_win{}'.format(pdir, nspec, win)
    if os.getcwd().startswith('/Users/davidmurphy'):
        pdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast'
        pth = '{}/aligned{}_win{}'.format(pdir, nspec, win)

    # get all the count files for the chromosome into a list
    arrs = []
    for f in os.listdir(pth):
        if f.startswith(ch + '.') and f.endswith('subcnt.txt'):
            f_name = '{}/{}'.format(pth, f)
            arrs.append(np.loadtxt(f_name))

    # convert list of counts to single array
    arrs = np.concatenate(arrs)
    # get unique sorting index from start positions
    sidx = np.unique(arrs[:,0], return_index=True)[1]
    arrs = arrs[sidx]

    # update first and last positions to cover the chromosome
    arrs[0, 0] = 0
    arrs[-1, 1] = chromosome_length(ch)

    # walk through the segments and fill any gaps
    f_save = '{}/{}.subcount.filled.txt'.format(pth, ch)
    rfmt = '{:.0f} {:.0f} {:.0f} {:.0f} {:.0f}  {:.0f}\n'
    n_missing = 0
    with open(f_save, 'w') as f:
        for i in range(len(arrs)-1):
            # get current line
            start1, end1, bases1, subs1, cpg1, bgc1 = arrs[i]
            # get following line
            start2, end2, bases2, subs2, cpg2, bgc2 = arrs[i+1]

            # if following line is noncontiguous, fill
            if start2 > end1:
                startx, endx = end1, start2
                # use average for interpolating count
                subx = int((subs1+subs2) / 2.0)
                cpgx = int((cpg1 + cpg2) / 2.0)
                bgcx = int((bgc1 + bgc2) / 2.0)
                # write current line
                f.write(rfmt.format(start1, end1, bases1, subs1, cpg1, bgc1))
                # write fill line
                f.write(rfmt.format(startx, endx, 7000, subx, cpgx, bgcx))
                n_missing += 1
            # if following line is contiguous, record and continue
            else:
                f.write(rfmt.format(start1, end1, bases1, subs1, cpg1, bgc1))

        # add final line
        f.write(rfmt.format(*arrs[-1]))
        # print the number of missing lines filled
        msg = 'filled missing lines = {}\n'.format(n_missing)
        stderr.write(msg)
        stdout.flush()


def main():
    spc_8 = ['hg19', 'gorGor3', 'ponAbe2', 'panTro4', 'papHam1', 'chlSab1',
             'macFas5', 'rheMac3']
    if len(argv) == 2:
        # ch = argv[1]
        # create_neutral_alignment(ch, spc_8)
        bbin = argv[1]
        # fix_chr20()
        # create_bbin_alignment(bbin)
        # process_bbin_subalign(bbin)
        # bbin_subcounts()
        # check_bbin()
    elif len(argv) == 3:
        ch, win = argv[1:]
        win = int(win)
        combine_subcount_files(ch, win)
    elif len(argv) == 5:
        ch = argv[1]
        start, end, win = map(int, map(float, argv[2:]))
        process_neutral_segment(ch, start, end, win)
    else:
        print('usage_1: estimate_substitution_rates <ch>')
        print('usage_2: estimate_substitution_rates <ch> <win>')
        print('usage_3: estimate_substitution_rates <ch> <start> <end> <win>')
        exit(1)


def main_2():
    if len(argv) != 3:
        print('usage: estimate_substitution_rates <folder> <bbin>')
    fldr, bbin = argv[1:]
    cst = cst_from_fldr(fldr)
    create_bbin_alignment(cst, bbin)


if __name__ == '__main__':
    main_2()
    # for c in range(2, 23):
    #     ch = 'chr{}'.format(c)
    #     combine_subcount_files(ch, 7000)
