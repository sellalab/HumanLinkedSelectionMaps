__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv, stderr, stdout
from classes.geneticmap import GeneticMap
from classes.runstruct import root_dir, cst_from_fldr
from data_processing.data_tools import wigz_array
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import load_saved, adjust_arrays


final_dir = root_dir + '/result/final_files'


def err_msg(errstring):
    stderr.write(errstring + '\n')
    stdout.flush()


def get_chrom_pred(cst):
    """return predictions for a single chromosome"""
    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst, chroms=[cst.chrom])

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    cutoff = 0.25
    mnu1 = (nu <= cutoff)

    # mask and rescale maps
    bs, cs, nu, nt, dv, pl, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    mnu2 = (nu <= cutoff)
    if bs is not None:
        bs = bs[mnu2]
    if cs is not None:
        cs = cs[mnu2]
    nu, nt, dv, pl = [a[mnu2] for a in nu, nt, dv, pl]
    msk &= mnu1[:, 0]

    # set uniform mutation rate
    nu_pred = np.full(len(nu), cst.stat.meandiv)

    # get predicted pi
    pred = predicted_pi(cst.params, cst, nu_pred, bs, cs)

    # get cumulative positions and mask them
    pos = np.cumsum(sg)[msk]
    assert len(pos) == len(pred)

    return pos, pred


def get_chrom_cons(cst):
    """get the phastcons scores for the neutral sites in nmsk"""
    # set fil path
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    fstr = '{p}/{cn}/{cn}.{c}.wig.gz'
    fwigz = fstr.format(p=pth, cn='fish', c=cst.chrom)

    # get conservation scores for fish, use -1 for missing
    con = wigz_array(cst.chrom, fwigz, use_neg1=True)

    return con


def indexed_positions(pos, idx, l=1e6):
    """get the region of the positions to use based on the index"""
    # define start index and end index for a given length
    i = int(idx*l)
    j = int(i+l)
    # if the index is too large for the neutral positions, exit with error msg
    if i >= len(pos):
        err_msg('index exceeds length of neutral positions.')
        exit(1)
    # otherwise return the subset of positions for the given index
    else:
        return pos[i:j]


def predsort_cons(ch, fldr, idx):
    """looking at conservation scores in vicinity of neutral sites on chrom"""

    # initialize ChromStruct
    cst = cst_from_fldr(fldr, ch)
    err_msg('cst loaded')

    # load neutral site positions for the chrom
    pos = np.where(np.load(cst.neut_masks)['neutmask'] > 0)[0]
    err_msg('pos loaded. len={}'.format(len(pos)))
    pos = indexed_positions(pos, idx)
    err_msg('indexed pos set. len={}'.format(len(pos)))

    # load predicted diversity levels for the chrom
    ipos, prd = get_chrom_pred(cst)
    err_msg('ipos and pred loaded')

    # load genetic map for the chrom
    gmp = GeneticMap(ch, cst.gmap_files)
    err_msg('gmap loaded')

    # load array of 99-vert conservation scores for the chrom
    con = get_chrom_cons(cst)
    err_msg('cons loaded')

    # create an empty array to store predictions, site counts and mean cons
    arr = np.zeros(shape=(len(pos), 5))

    # measure progress with percentiles of total positions
    pct = len(pos) / 100

    # # set path to save file
    # f_save = root_dir + '/compress/sort_cons/{}.{}.sort_cons.txt'.format(ch, fldr)
    # with open(f_save, 'w') as f:

    # loop through neutral positions and do the following:
    for i in xrange(len(pos)):
        # print out percentiles of progress
        if (i>0) and (i%pct == 0):
            if i/pct in [25, 50, 75]:
                stderr.write('{}%'.format(i/pct))
            else:
                stderr.write('.')
            stdout.flush()

        # get the current physical neutral site position
        p = pos[i]

        # interpolate the prediction at the site
        pred = np.interp(p, ipos, prd)

        # get genetic map position of neutral position
        gp = gmp.interp_gpos(p)

        # get genetic map position +/- 0.05cM from neutral position
        gp_i = max(0, gp-0.05)
        gp_j = gp+0.05

        # convert positions +/- 0.05cM from neutral position to physical pos
        p_i, p_j = gmp.interp_pos([gp_i, gp_j])

        # get all 99-vert conservation scores within this window
        cons_window = con[p_i:p_j]
        # get mask for sites with cons scores in the window
        wmask = (cons_window != -1)

        # store the predictions, number of sites and mean cons in window
        gmap_n_sites = np.sum(wmask)
        gmap_mean_con = np.mean(cons_window[wmask])

        # now get conservation scores in a 50 kb radius of the neutral site
        cons_window = con[p-50000:p+50000]
        # get mask for sites with scores in the window
        wmask = (cons_window != -1)
        kb_n_sites = np.sum(wmask)
        kb_mean_con = np.mean(cons_window[wmask])

        arr[i] = pred, gmap_n_sites, gmap_mean_con, kb_n_sites, kb_mean_con

        # # write output to file in real time
        # line = '{} {} {}\n'.format(pred, n_sites, mean_con)
        # f.write(line)

    # save the array of predictions, site counts and mean cons for the index
    fstr = root_dir + '/compress/sort_cons/{}.{}.{}.sort_cons.npz'
    f_save = fstr.format(ch, fldr, idx)
    np.savez_compressed(f_save, arr=arr)

    return None


def mean_in_bins(arr, nbins):
    """get the mean of the array values in bins"""
    step = int(len(arr) / nbins)
    mean_arr = []
    for i in xrange(0, len(arr), step):
        a = np.average(arr[i:i+step], axis=0)
        mean_arr.append(a)

    return np.array(mean_arr)


def combine_predsort_segments(fldr):
    """combine indexed sort data and average in windows"""
    # get all files across chromosomes and indices -- order does not matter
    sdir = root_dir + '/compress/sort_cons/'
    flist = []
    for f in os.listdir(sdir):
        if (fldr in f) and (f.endswith('sort_cons.npz')):
            flist.append(sdir + f)

    msg = 'found sort_cons {} files. preparing to load.'.format(len(flist))
    err_msg(msg)

    # iterate over files and add to list of arrays
    arr = []
    for f in flist:
        a = np.load(f)['arr']
        arr.append(a)
    err_msg('finished loading arrays into list')

    # concatenate arrays
    arr = np.concatenate(arr)
    err_msg('finished concatenating arrays')

    # get sorting index based on predictions and sort array
    arr = arr[np.argsort(arr[:,0])]
    err_msg('finished sorting arrays')

    # take mean in 100, 250, 500, 1000, 2000 windows
    for n in [100, 250, 500, 1000, 2000]:
        err_msg('taking mean values across {} bins'.format(n))
        mean_arr = mean_in_bins(arr, n)
        f_save = final_dir + '/{}/sorted.cons.n{}.txt'.format(fldr, n)
        np.savetxt(f_save, mean_arr)


def main():
    if len(argv) == 4:
        # load pre-saved data
        chrom, folder, idx = argv[1:]
        err_msg('chrom={} folder={} idx={} loaded.'.format(chrom, folder, idx))
        predsort_cons(chrom, folder, int(idx))
    elif len(argv) == 2:
        folder = argv[1]
        combine_predsort_segments(folder)
    else:
        print 'usage_1: prediction_sorted_conservation <chrom> <folder> <idx>'
        print 'usage_2: prediction_sorted_conservation <folder>'
        exit(1)

if __name__ == '__main__':
    main()
