#!/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python


__author__ = 'davidmurphy'


import os
from itertools import izip
from sys import argv, stderr, stdout
from classes.phylotree import parse_exptotsub
from classes.runstruct import chromosome_length, np, root_dir, ChromStruct


dloc = root_dir + '/data/phast/'
drem = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs/'


def err_msg(errstring):
    stderr.write(errstring + '\n')
    stdout.flush()


def outlier_rates(ch):
    """calculate substitution rates from exptotsub style files in a folder"""
    # chrom + "." for string match
    cl = ch + '.'
    if os.getcwd().startswith('/Users/davidmurphy'):
        din = dloc
    else:
        din = drem

    flist = [din+f for f in os.listdir(din) if cl in f and f.endswith('tsub')]

    # extra rates from every file and append to mast text file
    for fname in flist:
        # get chrom and base range by splitting filename
        ch, rng = fname.split('/')[-1].split('.')[:2]
        # get start, end positions from the range
        start, end = map(int, rng.split('_'))

        # process the expected substitution matrices from file
        dmat = parse_exptotsub(fname)
        # get the matrices for all human-macaque branches
        mats = [dmat[k] for k in dmat.keys() if ('hg19' in k or 'rheMac3' in k)]

        # calculate the total number of bases across branches
        bases = sum(m.sum() for m in mats)
        # calculate total substitutions across branches
        subs = bases - sum(m.diagonal().sum() for m in mats)

        # adjust bases by number of trees
        bases /= len(mats)

        # record filenames where subs >= bases
        if 0.5 * subs >= bases:
            print '{}.{}'.format(ch, rng), bases


def process_matrix(fname, return_string=True):
    """get sub counts from exptotsub matrix file and return formatted line"""
    # get chrom and base range by splitting filename
    ch, rng = fname.split('/')[-1].split('.')[:2]
    # get start, end positions from the range
    start, end = map(int, rng.split('_'))

    # process the expected substitution matrices from file
    dmat = parse_exptotsub(fname)
    # # get the matrices for all human-macaque branches
    # mats = [dmat[k] for k in dmat.keys() if ('hg19' in k or 'rheMac3' in k)]
    # get the matrices for ALL branches in the tree
    mats = dmat.values()

    # calculate the total number of bases across branches
    bases = sum(m.sum() for m in mats)
    # calculate total substitutions across branches
    subs = bases - sum(m.diagonal().sum() for m in mats)
    # calcuate CpG substitutions across branches
    cpg = sum(m[6, (4, 14)].sum() for m in mats)

    # adjust bases by number of trees
    try:
        bases /= len(mats)
    except ZeroDivisionError:
        err_msg('MATRIX FILE {} TRIGGERED ZERO DIVISION ERROR'.format(fname))

    # count total bases, subs and CpG adjusted subs
    b_cnt = int(bases)
    s_cnt = int(0.5 * subs)
    c_cnt = int(0.5 * (subs - cpg))

    # format results from current range to append to file
    fmt = '{st} {en} {bs} {sb} {cp}\n'
    line = fmt.format(st=start, en=end, bs=b_cnt, sb=s_cnt, cp=c_cnt)

    # return data string
    if return_string:
        return line
    # return tuple of data
    else:
        return start, end, b_cnt, s_cnt, c_cnt


def extract_rates(ch, fpath):
    """calculate substitution rates from exptotsub style files in a folder"""

    # determine the path based on current root
    if os.getcwd().startswith('/Users/davidmurphy'):
        din = dloc
    else:
        din = drem

    # set path to files
    pth = din + fpath + '/'
    # chrom + "." for string match
    cl = ch + '.'
    # get list of files matching the chromosome
    flist = [pth+f for f in os.listdir(pth) if cl in f and f.endswith('tsub')]

    # make path to output directory if it does not exist
    # dout = '/'.join(din.split('/')[:-2]) + '/sub_rates'
    dout = pth[:-1] + 'subcnt'
    if not os.path.isdir(dout):
        os.mkdir(dout)

    # # format output file name
    # fout = '{}/{}.sub_rates.txt'.format(dout, ch)
    # # add file header if file does not yet exist
    # if not os.path.isfile(fout):
    #     header = '#start end tot_rate noncpg_rate\n'
    #     with open(fout, 'w') as f:
    #         f.write(header)

    # format output file name
    fout = '{}/{}.sub_counts.txt'.format(dout, ch)
    # add file header if file does not yet exist
    if not os.path.isfile(fout):
        header = '#start end bases sub_count nonCpG_count\n'
        with open(fout, 'w') as f:
            f.write(header)

    # extra rates from every file and append to mast text file
    for fname in flist:
        # get chrom and base range by splitting filename
        ch, rng = fname.split('/')[-1].split('.')[:2]
        # get start, end positions from the range
        start, end = map(int, rng.split('_'))

        # process the expected substitution matrices from file
        dmat = parse_exptotsub(fname)
        # get the matrices for all human-macaque branches
        # mats = [dmat[k] for k in dmat.keys() if ('hg19' in k or 'rheMac3' in k)]
        # take all of the branches
        mats = dmat.values()

        # calculate the total number of bases across branches
        bases = sum(m.sum() for m in mats)
        # calculate total substitutions across branches
        subs = bases - sum(m.diagonal().sum() for m in mats)
        # calcuate CpG substitutions across branches
        cpg = sum(m[6,(4, 14)].sum() for m in mats)

        # adjust bases by number of trees
        bases /= len(mats)

        # # calculate total and CpG corrected substitution rates
        # # NOTE: multiply dimer rates by 0.5 to get per-bp (monomer) rates
        # tot_rate = 0.5 * subs / bases
        # cpg_rate = 0.5 * (subs-cpg) / bases
        #
        # # format the results for the current range (adjust start to 0-based)
        # fmt = '{st} {en} {tot} {cpg}\n'
        # line = fmt.format(st=start-1, en=end, tot=tot_rate, cpg=cpg_rate)

        # count total bases, subs and CpG adjusted subs
        b_cnt = int(bases)
        s_cnt = int(0.5 * subs)
        c_cnt = int(0.5 * (subs-cpg))

        # format results from current range to append to file
        fmt = '{st} {en} {bs} {sb} {cp}\n'
        line = fmt.format(st=start, en=end, bs=b_cnt, sb=s_cnt, cp=c_cnt)

        # append to the file containing all reults for the chrom
        with open(fout, 'a') as f:
            f.write(line)

    # # after appending all rates to the text file, save as npz with gaps filled
    # sdict = {}
    # # load sub count table into dict with start as keys
    # for (start, end, bases, subs, noncpg) in np.loadtxt(fout):
    #     sdict[int(start)] = [bases, subs, noncpg]
    #
    # # iterate over the chromosome in units of 20kb and fill gaps with nans
    # win = int(2e4)  # NOTE: use 20kb by default for now
    # new_table = []
    # for start in xrange(0, chromosome_length(ch), win):
    #     if start in sdict:
    #         row = [start] + sdict[start]
    #     else:
    #         row = [start, np.nan, np.nan, np.nan]
    #
    #     new_table.append(row)
    #
    # # convert table to array and save
    # new_table = np.array(new_table)
    # fsave = fout.replace('.txt', '.npz')
    # np.savez_compressed(fsave, rates=new_table)


def fill_gaps(ch):
    """fill in gaps with no data as nan for substitution rates"""
    # set path to subrate file
    din = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    fname = '{}/sub_rates/{}.sub_rates.txt'.format(din, ch)

    # load sub rate table into dict with start as keys
    sdict = {}
    for (start, end, rate, cpgrate) in np.loadtxt(fname):
        sdict[int(start)] = [rate, cpgrate]

    # iterate over the chromosome in units of 20kb and fill gaps with nans
    win = int(2e4)  # use 20kb by default for now
    new_table = []
    for start in xrange(0, chromosome_length(ch), win):
        if start in sdict:
            row = [start] + sdict[start]
        else:
            row = [start, np.nan, np.nan]

        new_table.append(row)

    # convert table to array and save
    new_table = np.array(new_table)
    fsave = fname.replace('.txt', '.npz')
    np.savez_compressed(fsave, rates=new_table)


def fill_gaps_2(fname, win):
    # get chrom from filename
    ch = fname.split('/')[-1].split('.')[0]
    # load sub counts for chromosome
    sarr = np.loadtxt(fname)
    # FIX FOR SLIDING WINDOWS BUG!! CUTOFF EXTRA WINDOW!!
    if np.isnan(sarr[-2, 3]) and np.isnan(sarr[-1, 3]):
        sarr = sarr[:-1]
        print 'fixed bad end'
    # filter out problematic last segments
    msk = sarr[:, 2] < 1.25 * win
    start, end, bases, sub, cpgsub = sarr[msk].T
    # convert start/end to int for indexing
    start, end = [a.astype(int) for a in start, end]
    # check that coordinates make sense
    assert np.all(end[:-1] == start[1:])

    # reset start and end to fill the chromosome
    start[0] = 0
    end[-1] = chromosome_length(ch)
    # check that counts now cover chromosome
    assert np.sum(end - start) == chromosome_length(ch)

    # fix nans
    m = np.isnan(sub)
    sub[m] = np.interp(start[m], start[~m], sub[~m]).astype(int)
    cpgsub[m] = np.interp(start[m], start[~m], cpgsub[~m]).astype(int)

    # rebuild and return array
    return np.column_stack((start, end, bases, sub, cpgsub))


def join_counts(ch, fpath, nsites=1e4):
    """join the separate sub counts across files for a chromosome"""
    # add "." to chrom string to avoid ambiguous matches (i.e., chr2/chr22)
    cl = ch + '.'
    # add final "/" to fpath if needed
    if not fpath.endswith('/'):
        fpath += '/'

    # combine all of the matching files into one array
    joined = []
    for f in os.listdir(fpath):
        if (cl in f) and (f.endswith('subcnt.txt')):
            joined.append(np.loadtxt(fpath+f))
    joined = np.concatenate(joined)

    # sort the sub count segments using start positions
    joined = joined[np.argsort(joined[:,0])]

    # extend the first and last window to the start and end of the chrom
    joined[0, 0] = 0
    joined[-1, 1] = chromosome_length(ch)

    # check for any remaining gaps between end to start of segments
    new_segs = []
    for (start, end) in izip(joined[1:, 0], joined[:-1, 1]):
        if start != end:
            # create new segment to fill gap. use NaNs for subcount values
            new_segment = [end, start, nsites, np.nan, np.nan]
            new_segs.append(new_segment)

    # combine the new segments with the original segments and re-sort
    if new_segs:
        joined = np.concatenate((joined, np.array(new_segs)))
        joined = joined[np.argsort(joined[:, 0])]

    # check that all gaps have been filled
    assert np.sum(joined[:,1] - joined[:,0]) == chromosome_length(ch)

    return joined


def expand_rates(ch, nspec, win):
    """convert compressed coords to spatial coords in subcount files"""
    fdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    # constuct file path
    fpath = fdir + '/aligned_{}_win_{}_subcounts'.format(nspec, win)
    # construct filename
    fname_1 = fpath + '/{}.subcount.txt'.format(ch)
    fname_2 = fpath + '/{}.sliding.subcount.txt'.format(ch)

    for fname in fname_1, fname_2:
        new_fname = fname.replace('.txt', '.filled.txt')
        new_arr = fill_gaps_2(fname, win).astype(int)
        np.savetxt(new_fname, new_arr, fmt='%d %d %d %d %d')


def create_missing_list(fpath, opath):
    """create a set of files containing commands for re-running missing files"""
    nfiles = 0
    arg_list = []
    for afile in os.listdir(fpath):
        fname = fpath + '/' + afile
        # get ch, nspec and window from file name
        ch, nspec, win = afile.split('.')[:3]
        with open(fname, 'r') as f:
            for line in f:
                line = line.split()
                # process missing lines to add to commands list
                if line[3] == 'MISSING':
                    argline = '{} {} {} {}'.format(ch, line[1], nspec, win)
                    arg_list.append(argline)
                # every 600 lines, write list to file and empty list
                if len(arg_list) == 600:
                    ofile = '{}/missing.{}.txt'.format(opath, nfiles)
                    with open(ofile, 'w') as o:
                        o.write('\n'.join(arg_list))
                    nfiles += 1
                    arg_list = []


def find_missing(ch, nspec, win):
    """iterate through the expected filenames for run and find missing ones"""
    fdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    # constuct file path
    fpath = fdir + '/aligned_{}_win_{}'.format(nspec, win)

    # use neutral mask to get expected length of alignment
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = np.load(cst.neut_masks)['neutmask']
    alen = nmsk.sum()

    # chr10.10000000_10010000.exptotsub
    # run through expected files and print missing filenames
    shift = win / 2
    for i in xrange(0, alen, shift):
        j = i+win
        fname = fpath + '/{}.{}_{}.exptotsub'.format(ch, i, j)
        if not os.path.isfile(fname):
            estr = '{} {} {} MISSING'.format(ch, i, j)
            # estr = '{} {} {} {}'.format(ch, i, nspec, win)
        else:
            estr = '{} {} {} FOUND'.format(ch, i, j)

        # print error message to log
        err_msg(estr)


def find_missing_2(ch, nspec, win):
    """iterate through the expected filenames for run and find missing ones"""
    fdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    # fdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast'
    # constuct file path
    fpath = fdir + '/aligned_{}_win_{}'.format(nspec, win)

    # use neutral mask to get expected length of alignment
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = np.load(cst.neut_masks)['neutmask']
    alen = nmsk.sum()
    npos = np.where(nmsk == 1)[0]

    # combine PHYSICAL start positions from all of the files across chroms
    ipos = []
    for f in os.listdir(fpath):
        if ch + '.' in f:
            # print f
            ipos.append(np.loadtxt(fpath+'/'+f)[:,0])
    ipos = np.concatenate(ipos).astype(int)
    ipos.sort()

    # create index of neutral alignment start positions
    nidx = np.arange(0, alen, win).astype(int)
    # get the expected PHYSICAL start positions for the chromosome
    exp_ipos = npos[nidx]

    # get positions of missing expected start positions
    imissing = ~np.in1d(exp_ipos, ipos)
    # convert to missing neutral alignment start positions
    missing_idx = nidx[imissing]

    # create list of commands to rerun missing segments
    lnum = 0
    fnum = 0
    fout = open(fpath + '/{}.missing.{}.txt'.format(ch, fnum), 'w')
    for i in missing_idx:
        msg = '{} {} {} {}\n'.format(ch, i, nspec, win)
        fout.write(msg)
        lnum += 1
        # create new file every 500 lines
        if lnum % 500 == 0:
            fout.close()
            fnum += 1
            fout = open(fpath + '/{}.missing.{}.txt'.format(ch, fnum), 'w')

    # close last opened file
    fout.close()


def find_missing_3(ch, nspec, win, nlines=500):
    """iterate through the expected filenames for run and find missing ones"""
    # constuct file path
    if os.getcwd().startswith('/Users/davidmurphy/'):
        fdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast'
    else:
        fdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    pth = fdir + '/aligned_{}_win_{}'.format(nspec, win)

    # load start positions from combined subcount file
    ipos = np.loadtxt('{}/{}.subcount.txt'.format(pth, ch))[:,0]

    # create index of neutral alignment start positions from neutral mask
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = np.load(cst.neut_masks)['neutmask']
    npos = np.where(nmsk == 1)[0]
    step = win / 10
    nidx = np.arange(0, npos.size, step).astype(int)

    # get the expected PHYSICAL start positions for the chromosome
    exp_ipos = npos[nidx]
    # get positions of missing expected start positions
    imissing = ~np.in1d(exp_ipos, ipos)
    # convert to missing neutral alignment start positions
    missing_idx = nidx[imissing]

    # create list of commands to rerun missing segments
    if len(missing_idx):
        lnum = 0
        fnum = 0
        fout = open(pth + '/{}.missing.{}.txt'.format(ch, fnum), 'w')
        for i in missing_idx:
            msg = '{} {} {} {}\n'.format(ch, i, nspec, win)
            fout.write(msg)
            lnum += 1
            # create new file every 500 lines
            if lnum % nlines == 0:
                fout.close()
                fnum += 1
                fout = open(pth + '/{}.missing.{}.txt'.format(ch, fnum), 'w')

        # close last opened file
        fout.close()


def interp_missing(ch, nspec, win, dstep=10):
    """interpolate values for final missing segments"""
    # constuct file path
    if os.getcwd().startswith('/Users/davidmurphy/'):
        fdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast'
    else:
        fdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    pth = fdir + '/aligned_{}_win_{}'.format(nspec, win)

    # load start positions from combined subcount file
    sarr = np.loadtxt('{}/{}.subcount.txt'.format(pth, ch))
    ipos = sarr[:,0]

    # create index of neutral alignment start positions from neutral mask
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = np.load(cst.neut_masks)['neutmask']
    npos = np.where(nmsk == 1)[0]
    step = win / dstep
    nidx = np.arange(0, npos.size, step).astype(int)

    # create inner window start/end positions
    inner_start = (win / 2) - (step / 2)
    i_in = np.minimum(nidx+inner_start, npos.size-1)
    j_in = np.minimum(i_in+step, npos.size-1)
    assert nidx.size == i_in.size
    # get the PHYSICAL start/end positions for inner windows
    pi_in = npos[i_in]
    pj_in = npos[j_in]
    # adjust first start and final end to span chrom
    pi_in[0] = 0
    pj_in[-1] = chromosome_length(ch)

    # get the expected PHYSICAL start/end positions for the chromosome
    exp_ipos = npos[nidx]
    njdx = np.minimum(len(npos)-1, nidx+win)
    exp_jpos = npos[njdx]
    # get positions of missing expected start positions
    imissing = ~np.in1d(exp_ipos, ipos)
    # convert to missing neutral alignment start positions
    missing_pos = exp_ipos[imissing]

    # interpolate missing counts
    s_cnt, c_cnt = sarr[:,3], sarr[:,4]
    s_itp = np.interp(missing_pos, ipos, s_cnt)
    c_itp = np.interp(missing_pos, ipos, c_cnt)

    # create lines for interpolated data
    i, j = missing_pos, exp_jpos[imissing]
    bases = np.full(len(i), win)
    rows = np.column_stack((i, j, bases, s_itp, c_itp))

    new = np.concatenate((sarr, rows))
    sidx = np.argsort(new[:,0])
    new = new[sidx]

    # now ONLY keep expected rows (throws out smaller window shifts)
    eidx = np.in1d(new[:, 0], exp_ipos)
    new = new[eidx]

    # # get the half distance from start to end for each segment
    # hdis = (new[:,1] - new[:,0]).astype(int) / 2
    # # get the midpoint for each sliding window
    # mid = new[:,0] + hdis
    # # create new set of non-overlapping start/end points from mid points
    # inew = np.concatenate(([0], mid[:-1]))
    # jnew = np.concatenate((mid[:-1], [chromosome_length(ch)]))

    # update start/end points in new array
    new[:,0] = pi_in
    new[:,1] = pj_in

    assert np.all(new[:-1, 1] == new[1:, 0])
    # # reset start and end to fill the chromosome
    # new[0, 0] = 0
    # new[-1,1] = chromosome_length(ch)
    # check that counts now cover chromosome
    assert np.sum(new[:,1] - new[:,0]) == chromosome_length(ch)

    pct = int(100 * 1.0 / dstep)
    if dstep == 1:
        fout = '{}/{}.subcount.filled.txt'.format(pth, ch, pct)
    else:
        fout = '{}/{}.subcount.filled.{}pct.txt'.format(pth, ch, pct)
    np.savetxt(fout, new, fmt='%d %d %d %d %d')


def build_subcount_files(ch, nspec, win):
    """iterate through the expected filenames for run and find missing ones"""
    # base dir for phyloFit runs
    basedir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    # constuct phyloFit run file path
    fpath = basedir + '/aligned_{}_win_{}'.format(nspec, win)
    # construct output file path, create dir if needed
    opath = fpath + '_subcounts'
    if not os.path.isdir(opath):
        os.mkdir(opath)

    # use neutral mask to get expected length of alignment
    token = 'basic.ncons.aligned.{}'.format(nspec)
    cst = ChromStruct(chrom=ch, tkn=token)
    nmsk = np.load(cst.neut_masks)['neutmask']
    alen = nmsk.sum()
    npos = np.where(nmsk != 0)[0]

    # create two output files for nonoverlapping and sliding windows
    fout_1 = '{}/{}.subcount.txt'.format(opath, ch)
    fout_2 = '{}/{}.sliding.subcount.txt'.format(opath, ch)

    # set the format for output rows to files
    fmt = '{} {} {} {} {}\n'

    # run through expected files and process them
    shift = win / 2
    for i in xrange(0, alen, shift):
        # if bottom coordinate is larger than neutral pos, break
        if i > len(npos):
            err_msg('ERROR: INDEX={} EXCEEDS NEUTRAL MATRIX'.format(i))
            break
        # j = i+win
        # if top coordinate is larger than neutral pos, use min
        j = min(len(npos-1), i+win)

        # construct filename for given coords and find
        matfile = '{}.{}_{}.exptotsub'.format(ch, i, j)
        fname = fpath + '/' + matfile
        # if the expected file exists, process exptotsub matrix
        if os.path.isfile(fname):
            start, end, b_cnt, s_cnt, c_cnt = process_matrix(fname)
            assert (start == i) and (end == j)
        # otherwise create a filler line with nans
        else:
            err_msg('MISSING MATRIX FILE: {}'.format(matfile))
            start, end, b_cnt, s_cnt, c_cnt = i, j, win, np.nan, np.nan

        # for sliding windows, process every line and adjust values
        islide, jslide = i + shift/2, j - shift/2
        # convert indices to physical positions
        ipos, jpos = npos[islide], npos[jslide]
        # write values to sliding window file
        with open(fout_2, 'a') as f:
            f.write(fmt.format(ipos, jpos, b_cnt, s_cnt, c_cnt))

        # for nonoverlapping windows, use nonoverlapping values
        if i%win == 0:
            ipos, jpos = npos[i], npos[j]
            # write values to nonoverlapping window file
            with open(fout_1, 'a') as f:
                f.write(fmt.format(ipos, jpos, b_cnt, s_cnt, c_cnt))


def comb_subcount_files(ch, nspec, win):
    """combine original and "missing" files into single subcount"""
    # constuct file path
    if os.getcwd().startswith('/Users/davidmurphy/'):
        fdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast'
    else:
        fdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs'
    pth = fdir + '/aligned_{}_win_{}'.format(nspec, win)

    # check for existing combined subcount file
    fout = '{}/{}.subcount.txt'.format(pth, ch)
    if not os.path.isfile(fout):
        # first combine originals
        orig = []
        for f in os.listdir(pth):
            if (ch + '.' in f) and (f.endswith('subcnt.txt')):
                fpth = '{}/{}'.format(pth, f)
                orig.append(np.loadtxt(fpth))
                # os.remove(fpth)
        orig = np.concatenate(orig).astype(int)
    else:
        orig = np.loadtxt(fout)

    # next combine leftover 'missing' files
    miss = []
    for f in os.listdir(pth):
        fpth = '{}/{}'.format(pth, f)
        if f.startswith('{}.missing.'.format(ch)):
            if f.endswith('.rerun.txt'):
                miss.append(np.loadtxt(fpth))
            os.remove(fpth)

    # if there were any missing files, continue processing missing lines
    if len(miss):
        miss = np.concatenate(miss).astype(int)
        # convert to PHYSICAL positions in missing files using neutral mask
        token = 'basic.ncons.aligned.{}'.format(nspec)
        cst = ChromStruct(chrom=ch, tkn=token)
        nmsk = np.load(cst.neut_masks)['neutmask']
        npos = np.where(nmsk == 1)[0]
        # reset positions in the array
        miss[:,0] = npos[miss[:,0]]
        miss[:,1] = npos[miss[:,1]]
        # combine original and missing lines into single array and sort
        comb = np.concatenate((orig, miss))
    else:
        comb = orig

    sidx = np.argsort(comb[:,0])
    comb = comb[sidx]

    # save combined data into single output file
    fout = '{}/{}.subcount.txt'.format(pth, ch)
    np.savetxt(fout, comb, fmt='%d %d %d %d %d')


if os.getcwd().startswith('/Users/davidmurphy'):
    for ds in [1]:
        for c in xrange(1, 23):
            ch = 'chr{}'.format(c)
            interp_missing(ch, 8, 7000, ds)
            print '{} {}% shift done.'.format(ch, 100 / ds)


def main():
    if len(argv) != 4:
        print 'usage: get_subrates <chrom> <nspec> <win>'
        exit(1)

    # get all the files from the current chromosome
    # fill_gaps(ch)
    # outlier_rates(ch)
    # extract_rates(ch, fpath)
    ch = argv[1]
    nspec, win = map(int, argv[2:])
    # find_missing_2(ch, nspec, win)
    # build_subcount_files(ch, nspec, win)
    # expand_rates(ch, nspec, win)
    # comb_subcount_files(ch, nspec, win)
    # find_missing_3(ch, nspec, win, nlines=50)
    interp_missing(ch, nspec, win)


if __name__ == '__main__':
    if os.getcwd().startswith('/ifs/'):
        main()
