__author__ = 'davidmurphy'


import os
import gzip
import numpy as np
from sys import stderr, stdout
from neutral_data import FilePaths, get_neutral_masks
from classes.runstruct import ChromStruct


cm_bins = np.concatenate(([1e-4, 1e-3], np.arange(0.01, 1.01, 0.01), [10]))

"""
BOOTSTRAP TOOLS
"""


def window_indices(token, aut=False, wsize=1e6):
    """get indices that divide X/A data into <wsize> windows"""
    # initialize a chromstruct
    cst = ChromStruct('chr1', outg='rheMac3', tkn=token)
    fpth = FilePaths(token)

    # dividing chrX is trivial
    if not aut:
        # load the neutral mask array: non-zero positions are neutral sites
        cst.chrom = 'chrX'
        nmsk = get_neutral_masks(cst)

        npos = np.where(nmsk != 0)[0] + 1  # convert to 1-based coordinates

        # sanity check that positions match arrays used for X/A analyses
        assert npos.size == np.load(fpth.fbval('chrX'))['bval'].size

        # adjust final window endpoint to cover all bases in chrX
        end = cst.chlen + wsize - cst.chlen % wsize + 1

    # for autosomes, we divide concatenated chromosomes into 1MB bins
    else:
        aut_npos = []
        aut_len = 0
        for ch in cst.chroms:
            # get neutral positions from neutral mask
            cst.chrom = ch
            nmsk = get_neutral_masks(cst)
            npos = np.where(nmsk != 0)[0] + 1  # convert to 1-based coordinates

            # check against file used in X/A analyses
            assert npos.size == np.load(fpth.fbval(ch))['bval'].size

            # augment npos by cumulative sum of previous chrom lengths
            npos += aut_len

            # increment the total cumulative sum by current chrom
            aut_len += cst.chlen

            # record neutral positions
            aut_npos.append(npos.astype('uint32'))

        # join the cumulative autosomal positions
        npos = np.concatenate(aut_npos)

        # adjust final window endpoint to cover all bases in joined autosomes
        end = aut_len + wsize - aut_len % wsize + 1

    # find indices that sort neutral sites into 1MB intervals in X and Aut
    wins = np.arange(0, end, wsize, dtype='uint32')
    widx = np.searchsorted(npos, wins)

    return widx


def create_blocks(neut, token, outg='rheMac3', sortby='bval', aut=True, wsize=1e6):
    fpth = FilePaths(token)
    # NOTE: fixing to McVicker maps only
    mcv = True

    # set files for either X or A data
    if aut:
        bfile = fpth.abval_file if mcv else fpth.anewbval_file
        pifile = fpth.fapi(neut)
        dvfile = fpth.fadiv(outg)
        cmfile = fpth.adist_file
    else:
        bfile = fpth.fbval('chrX')
        pifile = fpth.fpi('chrX', neut)
        dvfile = fpth.fdiv('chrX', outg)
        cmfile = fpth.fdist('chrX')

    # load bvals for each neutral site
    bvals = np.load(bfile)['bval']

    if sortby == 'dist':
        # load distance to nearest gene in cm and sort into bins
        dists = np.load(cmfile)['dist']
        bin_num = np.searchsorted(cm_bins, dists)
        nbins = len(cm_bins)
        del dists
    else:
        bin_num = bvals
        nbins = 1001

    # load pi and div arrays
    api = np.load(pifile)['pi']
    adiv = np.load(dvfile)['div']

    # generate window indices
    widx = window_indices(token, aut, wsize)
    # bin data in "wsize" length blocks
    blk_list = []
    for (i, j) in zip(widx[:-1], widx[1:]):
        # get the pi values for the window
        pi_block = api[i:j]
        # get the d values for the window
        d_block = adiv[i:j]
        # get the b values for the window
        b_block = bvals[i:j]

        # get the bin labels for the window
        if sortby == 'dist':
            block = bin_num[i:j]
        else:
            block = b_block

        # create empty arrays for each window with bins for each b val
        cn = np.zeros(shape=nbins)
        pi = np.zeros(shape=nbins)
        dv = np.zeros(shape=nbins)
        bv = np.zeros(shape=nbins)

        # for each window, summarize counts, pi and div for each b val bin
        for ii in xrange(nbins):
            iidx = np.where(block == ii)[0]
            cn[ii] = len(iidx)
            pi[ii] = np.sum(pi_block[iidx])
            dv[ii] = np.sum(d_block[iidx])
            bv[ii] = np.sum(b_block[iidx])

        # combine the summary vectors into an array
        if sortby == 'dist':
            blk = np.column_stack((cn, pi, dv, bv))
        else:
            blk = np.column_stack((cn, pi, dv))

        blk_list.append(blk)

    # save all blocks as a single file
    if aut:
        f = fpth.fblck_list(neut, sortby, wsize)
    else:
        f = fpth.fxblck_list(neut, sortby, wsize)
    np.savez_compressed(f, blist=np.array(blk_list))

    return None


def combine_bins(blk, bvals, min=10000):
    # create a mask to remove b value bins with no data
    msk = blk[:, 0] > 0
    blk = blk[msk]
    bvals = bvals[msk]

    # merge rows until each row has >= min site counts
    low_cnt = (blk[:, 0] < min)
    while np.any(low_cnt):
        # get the index of the FIRST low data bin
        ilow = np.where(low_cnt)[0][0]
        # add the data to the next bin unless ilow falls in the last bin
        iadd = ilow+1 if ilow+1 < len(blk) else ilow-1

        # take weighted sum of b values for combined bins
        n1, n2 = blk[iadd,0], blk[ilow,0]
        new_b = (n1 * bvals[iadd] +  n2 * bvals[ilow]) / (n1+n2)

        # combine data from ilow to iadd bins
        blk[iadd] += blk[ilow]
        bvals[iadd] = new_b

        # create new array with low bins removed
        msk = np.ones(shape=len(blk), dtype=bool)
        msk[ilow] = False
        blk = blk[msk]
        bvals = bvals[msk]

        # re-check counts
        low_cnt = (blk[:, 0] < min)

    return bvals, blk


def resample_bsdata(neut, token, smpl=10000, bmin=600, bmax=900, aut=True, wsize=1e6):
    # get the block list file
    fpth = FilePaths(token)
    sortby = 'bval'
    if aut:
        fblocks = fpth.fblck_list(neut, sortby, wsize)
    else:
        fblocks = fpth.fxblck_list(neut, sortby, wsize)

    # load blocks list file
    blks = np.load(fblocks)['blist']

    # take only b values in range of [bmin, bmax] from each block
    blks = blks[:, bmin:bmax+1, :]

    # remove blocks where there is no data for any b value
    blks = blks[(np.sum(blks[:, :, 0], axis=1) != 0)]

    # create a vector of 0-1 scaled b values for the b range used
    bvals = np.arange(bmin, bmax+1, 1) / 1e3

    # create a list of 10K resample indices to resample the blocks
    n = len(blks)
    sidx = [np.random.randint(0, n, size=n) for _ in xrange(smpl)]

    # use indices to resample blocks and record overall pi/BD per resample
    resample_means = []
    stderr.write('summarizing bootstrapped data: ')
    stdout.flush()
    for (cnt, si) in enumerate(sidx):
        # sum site counts, pi and div across sampled blocks
        sblk = np.sum(blks[si], axis=0)

        # remove empty b value bins and combine sparsely populated bins
        sbvals, sblk = combine_bins(sblk, bvals)

        # split count, pi and div columns into separate variables
        scn, spi, sdv = sblk.T

        # calculate pi scaled by divergence*B for each b value bin
        pi_bd = spi / (sdv*sbvals)

        # calculate pi scaled by B for each b value bin
        pi_b = spi / (sbvals*scn)

        # record site count weighted average for pi/BD genome-wide
        mean_pi_bd = np.average(pi_bd, weights=scn)

        # record site count weighted average for pi/B
        mean_pi_b = np.average(pi_b, weights=scn)

        # record individual weighted components of pi_bd
        mean_pi = np.average(spi/scn, weights=scn)
        mean_d = np.average(sdv/scn, weights=scn)
        mean_b = np.average(sbvals, weights=scn)

        resample_means.append((mean_pi_bd, mean_pi_b, mean_pi, mean_d, mean_b))

        if cnt % 100 == 0:
            msg = '.' if cnt%(100 * 25) else '{}%'.format(cnt/100)
            stderr.write(msg)
            stdout.flush()

    stderr.write('\n')
    stdout.flush()

    # save the resampled mean pi/BD values as a text file
    f = fpth.fpibd(neut, bmin, bmax) if aut else fpth.fxpibd(neut, bmin, bmax)
    np.savez_compressed(f,  bstrp=np.array(resample_means))

    return None


def summarize_range(neut, token, bmin=600, bmax=900, aut=True, wsize=1e6):
    """use resample blocks to summarize pi/B for a range of B values"""
    # get the block list file
    fpth = FilePaths(token)
    sortby = 'bval'
    if aut:
        xaut = 'A'
        fblocks = fpth.fblck_list(neut, sortby, wsize)
    else:
        xaut = 'X'
        fblocks = fpth.fxblck_list(neut, sortby, wsize)

    # load blocks list file
    blks = np.load(fblocks)['blist']

    # take only b values in range of [bmin, bmax] from each block
    blks = blks[:, bmin:bmax+1, :]

    # remove blocks where there is no data for any b value
    blks = blks[(np.sum(blks[:, :, 0], axis=1) != 0)]

    # create a vector of 0-1 scaled b values for the b range used
    bvals = np.arange(bmin, bmax+1, 1) / 1e3

    # sum site counts, pi and div across all blocks
    sblk = np.sum(blks, axis=0)

    # remove empty b value bins and combine sparsely populated bins
    sbvals, sblk = combine_bins(sblk, bvals)

    # split count, pi and div columns into separate variables
    scn, spi, sdv = sblk.T

    # calculate pi scaled by divergence*B for each b value bin
    pi_bd = spi / (sdv * sbvals)

    # record site count weighted average for pi/BD genome-wide
    mean_pi_bd = np.average(pi_bd, weights=scn)

    # record individual weighted components of pi_bd
    mean_pi = np.average(spi / scn, weights=scn)
    mean_d = np.average(sdv / scn, weights=scn)
    mean_b = np.average(sbvals, weights=scn)

    mean_vals = mean_pi_bd, mean_pi, mean_d, mean_b

    # save mean values to text file
    f_save = fpth.f_allblk(neut, xaut, bmin, bmax)
    f_head = '#mean_pi_bd mean_pi mean_d mean_b\n'
    f_str = '{} {} {} {}\n'.format(*mean_vals)
    with open(f_save, 'w') as f:
        f.write(f_head + f_str)

    return None


def resample_cmdata(neut, token, nsamples=10000, imin=0, aut=True, wsize=1e6):
    # get blocks list file
    fpth = FilePaths(token)
    sortby = 'dist'
    if aut:
        fblocks = fpth.fblck_list(neut, sortby, wsize)
    else:
        fblocks = fpth.fxblck_list(neut, sortby, wsize)

    # load blocks list file
    blks = np.load(fblocks)['blist']

    # remove blocks where there is no data for any cm bin
    blks = blks[(np.sum(blks[:, :, 0], axis=1) != 0)]

    imax = len(cm_bins)
    if imin >= imax:
        msg = 'ERROR: min cm bin {} >= total cm bins {}'
        raise ValueError(msg)

    # use imin to remove sites from blocks below a cM threshold
    blks = blks[:, imin:imax, :]

    # generate random indices to resample blocks for each bootstrap
    n = len(blks)
    sidx = [np.random.randint(0, n, size=n) for _ in xrange(nsamples)]

    # use indices to resample blocks and record overall pi/BD per resample
    resample_means = []
    stderr.write('summarizing bootstrapped data: ')
    stdout.flush()
    for (cnt, si) in enumerate(sidx):
        # create a new genome of resampled blocks using random block indices
        sblk = blks[si]

        # sum counts for each bvalue in blocks
        scn = np.sum(sblk[:, :, 0], axis=0)

        # mask out b value bins with no data
        nd = (scn != 0)
        scn = scn[nd]

        # sum pi for each cm bin and divide by block count
        spi = np.sum(sblk[:, :, 1], axis=0)[nd] / scn
        # sum div for each cm bin and divide by block count
        sdv = np.sum(sblk[:, :, 2], axis=0)[nd] / scn
        # sum bvals for each cm bin and divide by block count
        sbv = np.sum(sblk[:, :, 3], axis=0)[nd] / scn

        if (sdv == 0).any():
            print 'WARNING: D=0 encountered.'

        # calculate pi scaled by divergence for each cM bin
        pi_d = spi / sdv

        # record the weighted average value for pi/D for the whole sample
        mean_pi_d = np.average(pi_d, weights=scn)
        mean_pi = np.average(spi, weights=scn)
        mean_dv = np.average(sdv, weights=scn)
        mean_b = np.average(sbv, weights=scn)

        res = (mean_pi_d, mean_pi, mean_dv, mean_b, scn.sum())

        resample_means.append(res)

        if cnt % 100 == 0:
            msg = '.' if cnt % (100 * 25) else '{}%'.format(cnt / 100)
            stderr.write(msg)
            stdout.flush()

    stderr.write('\n')
    stdout.flush()

    # create a threshold label: cM in .3e format
    thresh_bins = np.concatenate(([0], cm_bins))
    thresh = '{:.3e}'.format(thresh_bins[imin])

    # save array of resampled mean pi/D values as a text file
    if aut:
        fsave = fpth.fpid(neut, thresh, wsize)
    else:
        fsave = fpth.fxpid(neut, thresh, wsize)

    # create a header
    hdr = 'pi/D\tpi\tD\tB\tsites'
    np.savetxt(fsave, np.array(resample_means), header=hdr, delimiter='\t')

    return None


def main():
    # run locally to debug
    if os.getcwd().startswith('/Users/davidmurphy'):
        # new_aut_bmap()
        # 100kb, 500kb, 2Mb, 5Mb
        neut = 'YRI'
        token = 'xaut_rmbgc'
        sortby = 'dist'
        aut = 1
        wsize = 1e6
        # resample_cmdata(neut, imin=32, aut=aut, wsize=wsize)
        create_blocks(neut, token, sortby='dist')
    # run on cluster
    else:
        from sys import argv
        # create blocks for X/A
        if len(argv) == 5:
            neut, token, sortby = argv[1:4]
            wsize = float(argv[4])
            create_blocks(neut, token, sortby=sortby, aut=0, wsize=wsize)
            create_blocks(neut, token, sortby=sortby, aut=1, wsize=wsize)
        # resample cM blocks
        elif len(argv) == 6:
            neut = argv[1]
            token = argv[2]
            imin = int(argv[3])
            aut = eval(argv[4])
            wsize = float(argv[5])
            resample_cmdata(neut, token, imin=imin, aut=aut, wsize=wsize)
        # resample B blocks
        elif len(argv) == 7:
            neut = argv[1]
            token = argv[2]
            bmin, bmax = int(argv[3]), int(argv[4])
            aut = eval(argv[5])
            wsize = float(argv[6])
            resample_bsdata(neut, token, bmin=bmin, bmax=bmax, aut=aut,
                            wsize=wsize)
        else:
            print 'usage_1: create_blocks <neut> <sortby> <wsize>'
            print 'usage_2: resample_cmdata <neut> <imin> <aut> <wsize>'


if __name__ == '__main__':
    main()
    # main1()
