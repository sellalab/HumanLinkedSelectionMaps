__author__ = 'davidmurphy'


import os
import subprocess
import numpy as np
from sys import argv, stderr, stdout
from datetime import datetime as dtime
from get_subrates import process_matrix
import data_processing.data_tools as dtl
from classes.runstruct import ChromStruct


def err_msg(errstring):
    stderr.write(errstring + '\n')
    stdout.flush()


def build_neutmask(cst):
    """get neutral mask based on conservation, sequence calling and filters"""

    # LOAD DATA
    # phastcons 100-vert p=0 sites
    ncmask = dtl.conservation_mask(cst.wigz_ncons, cst.chrom)
    # strict calling mask from 1000 genomes project
    callmask = dtl.sequence_mask(cst.call_mask, cst.chrom)
    # aggregate filter regions into a binary mask array
    filtmask = dtl.get_filter_mask(cst.chrom, cst.files.fnff)

    # BUILD MASK
    # initialize each position to True (i.e. neutral)
    neutmask = np.ones(shape=ncmask.size, dtype=bool)
    # default mask: 100-vertebrates phastCons score == 0
    neutmask &= ncmask
    # remove non-passing sites in polymorphism call mask (bool)
    neutmask &= callmask
    # apply optional filtering mask
    neutmask &= filtmask
    # save to filename
    np.savez_compressed(cst.neut_masks, neutmask=neutmask)

    return neutmask


def load_neutmask(cst):
    """create neutral mask or load from existing file"""
    if os.path.isfile(cst.neut_masks):
        return np.load(cst.neut_masks)['neutmask']
    else:
        err_msg('BUILDING NEUTRAL MASK')
        return build_neutmask(cst)


def save_features(features, fout):
    """save the start/end coords of neutral features for a window"""
    # save features in gff formatted file
    with open(fout, 'w') as f:
        fmt = 'hg19\tconservation\tprimate_nc\t{}\t{}\n'
        for (start, end) in features:
            f.write(fmt.format(start, end))


def count_nonbases(sequence):
    """count all of the non ACGT base characters in lines of FASTA data"""
    return sequence.count('*') + sequence.count('-') + sequence.count('N')


def verify_subalign(fname, features, frac):
    """assess the alignment quality at a set of features"""
    # calculate the number of feature bases in the sub-alignment
    feat_bases = sum(j-i for (i,j) in features)

    # dictionary for seqs and temp variables for processing file
    ma_dict = {}
    cur_lbl = ''
    cur_seq = []

    # read in the sub-alignment
    with open(fname, 'r') as f:
        for line in f:
            # process label lines
            if line.startswith('>'):
                # record the previous label and feature sequences
                if cur_lbl:
                    # join each line of sequence
                    seq = ''.join(cur_seq)
                    # slice features from the total sequence, count non-ACGT
                    nbc = sum(count_nonbases(seq[i:j]) for (i, j) in features)
                    ma_dict[cur_lbl] = nbc
                # reset current label and sequence
                cur_lbl = line[2:].strip('\n')
                cur_seq = []
            # process sequence lines
            else:
                cur_seq.append(line.strip('\n'))

    # analyze final sequence
    seq = ''.join(cur_seq)
    nbc = sum(count_nonbases(seq[i:j]) for (i, j) in features)
    ma_dict[cur_lbl] = nbc

    # return true if all species have at least 50% of expected feature bases
    assert ma_dict['hg19'] == 0
    nfrac = 1.0 - frac
    return all(v < nfrac*feat_bases for v in ma_dict.values())


def analyze_subalign(ch, start, nspec, win):
    # err_msg('starting analyze subalign')
    multi = 'primate'
    end = start+win
    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    cfg_dir = '{}/cfgs'.format(phast_dir)
    save_dir = '{}/runs/aligned_{}_win_{}'.format(phast_dir, nspec, win)

    # set input FASTA
    fa_in = '{}/{}.{}.neutral.aligned.{}.fa'.format(fa_dir, ch, multi, nspec)

    # set output sub-alignment FASTA file name
    fa_sub = fa_in.replace('.fa', '.{}_{}.fa'.format(start+1, end))
    # create W-kb sub-alignment using msa_view program
    fmt = 'msa_view {fin} --start {st} --end {en} --refidx 1 > {fout}'
    msa_cmd = fmt.format(fin=fa_in, st=start+1, en=end, fout=fa_sub)
    # call msa_view through the shell
    subprocess.call(msa_cmd, shell=True)
    # err_msg('finished calling msa_view')

    # set phyloFit command parameters
    tree = '{}/{}.{}.tree'.format(cfg_dir, multi, nspec)
    # use physical positions for save label preference
    pref = '{}/{}.{}_{}'.format(save_dir, ch, start, end)
    smod = 'U2S'
    flag = '-E -p MED -Z -q'
    # run phyloFit on the sub-alignment
    fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {lb} -s {sm} {fl}'
    phy_cmd = fmt.format(tr=tree, fa=fa_sub, lb=pref, sm=smod, fl=flag)
    # call phyloFit through the shell
    subprocess.call(phy_cmd, shell=True)
    # err_msg('finished calling phyloFit')

    # set phyloFit output file names
    mat_file = '{}.exptotsub'.format(pref)
    mod_file = '{}.mod'.format(pref)

    # process phyloFit file to get sub counts and append to out file
    scnt_line = process_matrix(mat_file)
    # err_msg('finished processing exptotsub file')

    # delete unused mod file
    os.remove(mod_file)
    # delete matrix file after processing
    os.remove(mat_file)
    # delete the sub-alignment
    os.remove(fa_sub)

    return scnt_line


def main_0():
    """
    this script calls tools from the PHAST suite to calculate context-specific
    substitution rates within a phylogenetic tree that has been split into
    windows. intermediate files are cleaned up and a final processed
    substitution counts file is created
    """
    if len(argv) != 6:
        print 'usage: split_subrates <ch> <multi> <start> <end> <window>'
        exit(1)

    # record the start time for runtime logging
    start_time = dtime.now()

    # get the chromosome, multiple alignment label and start/end coords, window
    ch, multi = argv[1], argv[2]
    start, end, win = map(int, map(float, argv[3:]))

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    ft_dir = '{}/features'.format(phast_dir)
    save_dir = '{}/runs/primate_windows_5'.format(phast_dir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # format output file name for substitution counts file
    fout = '{}/{}.{}_{}.subcnt.txt'.format(save_dir, ch, start, end)

    # set path to input FASTA file
    fa_in = '{}/{}.{}.fa'.format(fa_dir, ch, multi)

    # load neutral mask to slice window features
    cst = ChromStruct(chrom=ch, tkn='basic.ncons')
    nmsk = load_neutmask(cst)

    # for the start end range, create 20kb sub-alignments for phyloFit
    for i in range(start, end, win):
        # slice region of interest from mask (adjust for 0-based coords)
        j = i+win
        region = nmsk[i:j]

        # if there are too few neutral sites, continue to next region
        if np.sum(region) < 100:
            err_msg('--> SKIPPING REGION {}:{} - LOW DATA'.format(i, j))
            continue

        # set output sub-alignment FASTA file name
        fa_sub = '{}/{}.{}.{}_{}.fa'.format(fa_dir, ch, multi, i+1, j)
        # create W-kb sub-alignment using msa_view program
        fmt = 'msa_view {fin} --start {st} --end {en} --refidx 1 > {fout}'
        msa_cmd = fmt.format(fin=fa_in, st=i+1, en=j, fout=fa_sub)
        # call msa_view through the shell
        subprocess.call(msa_cmd, shell=True)

        # get neutral features from the region
        feats = dtl.binary_mask_segments(region)
        # if there is too little data across species, continue to next region
        frac = 0.75
        if not verify_subalign(fa_sub, feats, frac):
            err_msg('--> SKIPPING REGION {}:{} - POOR ALIGNMENT'.format(i, j))
            os.remove(fa_sub)
            continue

        # create neutral features file for sub-alignment using neutral region
        f_feat = '{}/{}.{}.{}_{}.gff'.format(ft_dir, ch, multi, i + 1, j)
        save_features(feats, f_feat)

        # set phyloFit command parameters
        tree = '{}/{}.tree'.format(cfg_dir, multi)
        pref = '{}/{}.{}_{}'.format(save_dir, ch, i+1, j)
        smod = 'U2S'
        flag = '-E -p MED -Z -C 1 -q'
        # run phyloFit on the sub-alignment
        fmt = 'phyloFit -t {} -i FASTA {} -o {} -s {} -g {} {}'
        phy_cmd = fmt.format(tree, fa_sub, pref, smod, f_feat, flag)
        # call phyloFit through the shell
        subprocess.call(phy_cmd, shell=True)

        # set phyloFit output file names
        mat_file = '{}.1.exptotsub'.format(pref)
        mod_file = '{}.primate_nc.mod'.format(pref)

        # process phyloFit file to get sub counts and append to out file
        scnt_line = process_matrix(mat_file)
        with open(fout, 'a') as f:
            f.write(scnt_line)

        # delete unused mod file
        os.remove(mod_file)
        # delete matrix file after processing
        os.remove(mat_file)
        # delete the sub-alignment
        os.remove(fa_sub)
        # delete gff features file
        os.remove(f_feat)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))


def main_1():
    """
    this script calls tools from the PHAST suite to calculate context-specific
    substitution rates within a phylogenetic tree that has been split into
    windows. intermediate files are cleaned up and a final processed
    substitution counts file is created
    """
    if len(argv) != 6:
        print 'usage: split_subrates <ch> <nspec> <start> <end> <window>'
        exit(1)

    # record the start time for runtime logging
    start_time = dtime.now()

    # fix "multi" to "primate"
    multi = 'primate'

    # get the chromosome, multiple alignment label and start/end coords, window
    ch, nspec = argv[1], argv[2]
    start, end, win = map(int, map(float, argv[3:]))

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    save_dir = '{}/runs/aligned_{}'.format(phast_dir, nspec)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # format output file name for substitution counts file
    fout = '{}/{}.{}_{}.subcnt.txt'.format(save_dir, ch, start, end)

    # set path to input FASTA file
    fa_in = '{}/{}.{}.neutral.aligned.{}.fa'.format(fa_dir, ch, multi, nspec)

    # load neutral mask to slice window features
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = load_neutmask(cst)

    # get positions of neutral sites from the neutral mask
    npos = np.where(nmsk == 1)[0]

    # for the start end range, create 20kb sub-alignments for phyloFit
    for i in range(start, end, win):
        # set upper bound for region; don't exceed total length of FASTA
        j = min(len(npos)-1, i+win)
        # get PHYSICAL positions of i & j for writing subcount file
        ipos, jpos = npos[i], npos[j]

        # set output sub-alignment FASTA file name
        fa_sub = fa_in.replace('.fa', '.{}_{}.fa'.format(i+1, j))
        # create W-kb sub-alignment using msa_view program
        fmt = 'msa_view {fin} --start {st} --end {en} --refidx 1 > {fout}'
        msa_cmd = fmt.format(fin=fa_in, st=i+1, en=j, fout=fa_sub)
        # call msa_view through the shell
        subprocess.call(msa_cmd, shell=True)

        # set phyloFit command parameters
        tree = '{}/{}.tree'.format(cfg_dir, multi)
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

        # process phyloFit file to get sub counts and append to out file
        scnt_line = process_matrix(mat_file)
        with open(fout, 'a') as f:
            f.write(scnt_line)

        # delete unused mod file
        os.remove(mod_file)
        # delete matrix file after processing
        os.remove(mat_file)
        # delete the sub-alignment
        os.remove(fa_sub)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))


def main_2():
    """
    this script calls tools from the PHAST suite to calculate context-specific
    substitution rates within a phylogenetic tree that has been split into
    windows. intermediate files are cleaned up and a final processed
    substitution counts file is created
    """
    if len(argv) != 6:
        print 'usage: split_subrates <ch> <nspec> <start> <end> <window>'
        exit(1)

    # record the start time for runtime logging
    start_time = dtime.now()

    # fix "multi" to "primate"
    multi = 'primate'

    # get the chromosome, multiple alignment label and start/end coords, window
    ch, nspec = argv[1], argv[2]
    start, end, win = map(int, map(float, argv[3:]))

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    save_dir = '{}/runs/aligned_{}_win_{}'.format(phast_dir, nspec, win)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # format output file name for substitution counts file
    fout = '{}/{}.{}_{}.subcnt.txt'.format(save_dir, ch, start, end)

    # set path to input FASTA file
    fa_in = '{}/{}.{}.neutral.aligned.{}.fa'.format(fa_dir, ch, multi, nspec)

    # load neutral mask to slice window features
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = load_neutmask(cst)

    # get positions of neutral sites from the neutral mask
    npos = np.where(nmsk == 1)[0]

    # for the start end range, create 20kb sub-alignments for phyloFit
    shift = win / 10
    for i in range(start, end, shift):
        # set upper bound for region; don't exceed total length of FASTA
        j = min(len(npos)-1, i+win)

        # set output sub-alignment FASTA file name
        fa_sub = fa_in.replace('.fa', '.{}_{}.fa'.format(i+1, j))
        # create W-kb sub-alignment using msa_view program
        fmt = 'msa_view {fin} --start {st} --end {en} --refidx 1 > {fout}'
        msa_cmd = fmt.format(fin=fa_in, st=i+1, en=j, fout=fa_sub)
        # call msa_view through the shell
        subprocess.call(msa_cmd, shell=True)

        # set phyloFit command parameters
        tree = '{}/{}.{}.tree'.format(cfg_dir, multi, nspec)
        # use PHYSICAL positions for save label preference
        pref = '{}/{}.{}_{}'.format(save_dir, ch, npos[i], npos[j])
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

        # process phyloFit file to get sub counts and append to out file
        scnt_line = process_matrix(mat_file)
        with open(fout, 'a') as f:
            f.write(scnt_line)

        # delete unused mod file
        os.remove(mod_file)
        # delete matrix file after processing
        os.remove(mat_file)
        # delete the sub-alignment
        os.remove(fa_sub)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))


def main_3():
    """
    this script calls tools from the PHAST suite to calculate context-specific
    substitution rates within a phylogenetic tree that has been split into
    windows. intermediate files are cleaned up and a final processed
    substitution counts file is created
    """
    if len(argv) != 6:
        print 'usage: split_subrates <ch> <nspec> <start> <end> <window>'
        exit(1)

    # record the start time for runtime logging
    start_time = dtime.now()

    # fix "multi" to "primate"
    multi = 'primate'

    # get the chromosome, multiple alignment label and start/end coords, window
    ch, nspec = argv[1], argv[2]
    start, end, win = map(int, map(float, argv[3:]))

    # set directories used in run
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    fa_dir = '{}/fa/{}'.format(phast_dir, multi)
    save_dir = '{}/runs/aligned_{}_win_{}'.format(phast_dir, nspec, win)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    cfg_dir = '{}/cfgs'.format(phast_dir)

    # # format output file name for substitution counts file
    # fout = '{}/{}.{}_{}.subcnt.txt'.format(save_dir, ch, start, end)

    # set path to input FASTA file
    fa_in = '{}/{}.{}.neutral.aligned.{}.fa'.format(fa_dir, ch, multi, nspec)

    # load neutral mask to slice window features
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = load_neutmask(cst)

    # get positions of neutral sites from the neutral mask
    npos = np.where(nmsk == 1)[0]

    # use steps that are 1/2 the window size for sliding windows
    step = win / 2
    # for the start end range, create N-kb sub-alignments for phyloFit
    for i in range(start, end, step):
        # set upper bound for region; don't exceed total length of FASTA
        j = min(len(npos)-1, i+win)

        # set output sub-alignment FASTA file name
        fa_sub = fa_in.replace('.fa', '.{}_{}.fa'.format(i+1, j))
        # create W-kb sub-alignment using msa_view program
        fmt = 'msa_view {fin} --start {st} --end {en} --refidx 1 > {fout}'
        msa_cmd = fmt.format(fin=fa_in, st=i+1, en=j, fout=fa_sub)
        # call msa_view through the shell
        subprocess.call(msa_cmd, shell=True)

        # set phyloFit command parameters
        tree = '{}/{}.{}.tree'.format(cfg_dir, multi, nspec)
        # use physical positions for save label preference
        pref = '{}/{}.{}_{}'.format(save_dir, ch, i, j)
        smod = 'U2S'
        flag = '-E -p MED -Z -q'
        # run phyloFit on the sub-alignment
        fmt = 'phyloFit -t {tr} -i FASTA {fa} -o {lb} -s {sm} {fl}'
        phy_cmd = fmt.format(tr=tree, fa=fa_sub, lb=pref, sm=smod, fl=flag)
        # call phyloFit through the shell
        subprocess.call(phy_cmd, shell=True)

        # # set phyloFit output file names
        # mat_file = '{}.exptotsub'.format(pref)
        # mod_file = '{}.mod'.format(pref)

        # # process phyloFit file to get sub counts and append to out file
        # scnt_line = process_matrix(mat_file)
        # with open(fout, 'a') as f:
        #     f.write(scnt_line)
        #
        # # delete unused mod file
        # os.remove(mod_file)
        # # delete matrix file after processing
        # os.remove(mat_file)
        # delete the sub-alignment
        os.remove(fa_sub)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))


def main_4():
    # record the start time for runtime logging
    start_time = dtime.now()
    if len(argv) != 2:
        print 'usage: split_subrates <list_file>'
        exit(1)

    # get the list file name from command line
    list_file = argv[1]
    # open output file using modified input file name
    out_file = list_file.replace('.txt', '.rerun.txt')

    # generate and execute commands from list of missing data
    with open(list_file, 'r') as f:
        for line in f:
            # err_msg(line.strip('\n'))
            line = line.split()
            ch = line[0]
            start, nspec, win = map(int, line[1:])
            # write output to file
            row = analyze_subalign(ch, start, nspec, win)
            with open(out_file, 'a') as fout:
                fout.write(row)

    # remove the "missing" file
    os.remove(list_file)

    # calculate the runtime and write to log
    run_time = dtime.now() - start_time
    err_msg('RUN TIME: {}'.format(run_time))


if __name__ == '__main__':
    main_4()
