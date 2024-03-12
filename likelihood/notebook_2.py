__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv
import data_processing.functions as fnc
import data_processing.data_tools as dtl
from classes.geneticmap import GeneticMap
from classes.annosegments import AnnoSegments
from simulations.neutpoly_sim import simulated_data
from classes.runstruct import ChromStruct, root_dir, izip
from likelihood.cllh_functions import serial_cllh, predicted_pi
from precalc.lh_inputs import load_saved, adjust_arrays, prepare_inputs


#%%
def cpg_overlaps(cst):
    # get C positions of all CpGs from ancestor
    icpg = dtl.cpg_mask(cst.chrom, cst.ancs_files)
    # get snp pos, ref, alt
    snpcnt = dtl.snpcount(cst.snp_files, returnbases=True)
    issnp, ref, alt = snpcnt[0], snpcnt[1], snpcnt[3]
    # create SNP pairs of ref/alt
    pairs = np.core.defchararray.add(ref, alt)
    # find CpG C positions overlapping C>T SNPs
    cti = np.in1d(icpg+1, issnp[pairs == 'CT'])
    # find CpG G positions overlapping G>A SNPs
    gai = np.in1d(icpg + 2, issnp[pairs == 'GA'])
    # count CpG sites and CpG SNPs
    ncpg, ncpgsnps = len(icpg), cti.sum() + gai.sum()
    pctsnp = 1.0 * ncpgsnps / ncpg
    print('{} {} {} {}'.format(cst.chrom, ncpg, ncpgsnps, pctsnp))


#%%
# print 'chrom CpG-sites CpG-SNPs'
# for c in range(1, 23):
#     ch = 'chr{}'.format(c)
#     cpg_overlaps(ChromStruct(ch))

#%%
def printouts_1():
    folders = 'datarun_000'.split()
    init_files = []
    print('running...')
    for fldr in folders:
        fdir = root_dir + '/result/final_files/{}/'.format(fldr)
        comp_file = [fdir+f for f in os.listdir(fdir) if 'composite' in f]
        assert len(comp_file) == 1
        fnc.swap_root(comp_file[0])
        init_files.append(comp_file[0])

    cstructs = [ChromStruct(chrom='chr1', init=ifile) for ifile in init_files]
    # cst = cstructs[0]
    for cst in cstructs:
        cst.vars.num_cores = 1
        params = cst.params
        ipt = prepare_inputs(cst)

        #%%
        print(serial_cllh(params, ipt))
        #
        # #%% load maps for chromstruct
        # cst = cstructs[3]
        # params = cst.params
        # sg, bs, cs, nu, nt, dv, pl = load_saved(cst)
        # bs, cs, nu, nt, dv, pl, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)
        #
        #%% get composite bmap with params
        bs = ipt.bs
        min_bsx = cst.fixed.min_bsx
        uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
        bwt = uvec / cst.fixed.u_fix
        bsx = np.exp(np.dot(bs, bwt))
        if min_bsx:
            print('min_bsx={}'.format(min_bsx))
            bsx = np.maximum(bsx, min_bsx)

        print(bsx.min())


def sortbin_edges(sites, numbins):
    """
    get the upper indices of sorted data array that divide
    data into bins of equal neutral site counts
    """

    # get number of sites needed so that numbins x numsites = total sites
    numsites = int(np.sum(sites) / numbins)

    # find indices that partition cumulative sorted site count into numsites
    cumsites = np.cumsum(sites)
    bounds = np.arange(numsites, cumsites[-1], numsites)

    # get the ending index for each partition
    jdx = list(np.searchsorted(a=cumsites, v=bounds))

    # return a list of (start, end) indices for each partition
    return zip([0] + jdx[:-1], jdx)


def calc_log_lh(hom, het, pred):
    return hom * np.log(1-pred) + het * np.log(pred)


def get_comp_file(folder_name):
    # find composite file in folder
    fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1
    f_init = fdir + flist[0]
    return f_init


def get_bsx(ipt, params=None):
    # cst = ChromStruct(chrom='chr1', init=init_file)
    # cst.vars.num_cores = 1
    cst = ipt.cst
    ipt = prepare_inputs(cst)
    bs = ipt.bs
    bsx = bsx_values(cst, bs, params)
    # min_bsx = cst.fixed.min_bsx
    # uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
    # bwt = uvec / cst.fixed.u_fix
    # bsx = np.exp(np.dot(bs, bwt))
    # if min_bsx:
    #     bsx = np.maximum(bsx, min_bsx)
    return bsx


def bsx_values(cst, bs, params=None):
    if params is None:
        params = cst.params
    min_bsx = cst.fixed.min_bsx
    uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
    bwt = uvec / cst.fixed.u_fix
    bsx = np.exp(np.dot(bs, bwt))
    if min_bsx:
        bsx = np.maximum(bsx, min_bsx)
    return bsx


def get_cllh(ipt, params=None):
    # cst = ChromStruct(chrom='chr1', init=init_file)
    # cst.vars.num_cores = 1
    # ipt = prepare_inputs(cst)
    if params is None:
        params = ipt.cst.params
    return serial_cllh(params, ipt)


def get_pred(ipt, params=None):
    # cst = ChromStruct(chrom='chr1', init=init_file)
    # cst.vars.num_cores = 1
    # ipt = prepare_inputs(cst)
    cst = ipt.cst
    if params is None:
        params = cst.params
    pii = predicted_pi(params, cst, ipt.u, ipt.bs, ipt.cs)
    return pii


def print_numsites(init_file):
    cst = ChromStruct(chrom='chr1', init=init_file)
    n = 0
    for ch in cst.chroms:
        cst.chrom = ch
        n += np.sum(np.load(cst.neut_masks)['neutmask'] > 0)
    print('{}: {}'.format(cst.tkn, n))


def observed_vs_expected_cllh(hom, het, bsx, pred, cnt, num):
    # get sorting indices based on prediction value
    si = np.argsort(bsx)

    # sort all data segments used for calculations
    hom, het, bsx, pred, cnt = hom[si], het[si], bsx[si], pred[si], cnt[si]

    # get start and end indices per partition
    idx = sortbin_edges(sites=cnt, numbins=num)

    # simulate data from predictions
    shom, shet = simulated_data(pred, cnt, 216).T

    # gather data summaries
    sorted_array = []
    for (i, j) in idx:
        # get data for current slice
        ihom, ihet, ishom, ishet = hom[i:j], het[i:j], shom[i:j], shet[i:j]
        ibsx, ipred, icnt = bsx[i:j], pred[i:j], cnt[i:j]

        # calculate mean bsx per bin
        mean_bsx = np.average(ibsx, weights=icnt)

        # calculate observed CLLH per bin with real data
        log_lh = np.sum(calc_log_lh(ihom, ihet, ipred))

        # calculate expected CLLH per bin with simulated data
        sim_log_lh = np.sum(calc_log_lh(ishom, ishet, ipred))

        # collect data into lista
        sorted_array.append([mean_bsx, log_lh, sim_log_lh])

    return np.array(sorted_array)


def obsd_vs_expt_tail_cons(cst, bth, params=None):
    ssize = 216
    assert isinstance(cst, ChromStruct)
    chrom = cst.chrom
    # load arrays for calculations
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst, [chrom])
    bs, cs, nu, _, _, _, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)
    # clear unused arrays
    del nt, dv, pl

    # use cumsum to get segment start/end points
    assert sg.sum() == cst.chlen
    sj = np.cumsum(sg, dtype=int)
    si = np.concatenate(([0], sj[:-1]))
    # reformat sg as Nx2 start/end segment array
    sg = np.column_stack((si, sj))
    # mask segments with no data
    sg = sg[msk]

    # generate prediction and bsx maps
    if params is None:
        pred = predicted_pi(cst.params, cst, nu, bs, cs)
    else:
        pred = predicted_pi(params, cst, nu, bs, cs)
    bsx = bsx_values(cst, bs, params)
    # mask data regions where bsx is <= threshold value
    bm = (bsx <= bth)
    pred, bsx, sg = pred[bm], bsx[bm], sg[bm]
    # clear unused arrays
    del nu, bs, cs

    # load neutral mask
    if os.getcwd().startswith('/Users/davidmurphy'):
        fnmsk = cst.neut_masks.replace('cleanrun', 'nff_2')
        b_anno = 'primate_cons95_Segments'
    else:
        fnmsk = cst.neut_masks
        b_anno = 'primate_cons94_Segments_rm_neut'
    print('USING {} AS DEFAULT BS ANNOTATIONS'.format(b_anno))
    nmsk = np.load(fnmsk)['neutmask']

    # load conserved segments as inverted mask (i.e., cons=0)
    cons = np.loadtxt(cst.bs_target(b_anno), usecols=(1, 2), dtype=int)
    cms = dtl.mask_segments(np.ones(shape=cst.chlen), cons, flipoff=True)

    # load SNP counts
    snpcnt = dtl.snpcount(cst.snp_files)

    # iterate over low B segments and find conserved regions
    segdata = []
    for (i, seg) in enumerate(sg):
        # get prediction and B values for current segment
        prd = np.array([pred[i]])
        # get the start and end of B segment
        si, sj = seg
        # create an array of positions within the segment
        spos = np.arange(si, sj, dtype=int)
        # slice the inverted cons mask at the current segment
        cons = cms[si:sj]
        # convert inverted cons mask to subsegments of current segment
        csegs = dtl.binary_mask_segments(cons)
        # slice the neutral mask at the current segment
        neut = nmsk[si:sj]

        # iterate over inverted cons segments, get neutral positions
        sdata = np.zeros(shape=(len(csegs), 5))
        for (ii, (csi, csj)) in enumerate(csegs):
            # add bsx column
            sdata[ii, 0] = bsx[i]
            # get the nmsk subsegment
            nslice = neut[csi:csj]
            # get neutral positions of the subsegment using nmsk slice vals
            pslice = spos[csi:csj]
            # simulate poly data for number of sites in subsegment
            ns = np.array([np.sum(nslice > 0)])
            sim = simulated_data(prd, ns, ssize, verbose=False)[0]
            sdata[ii, 1:3] = sim
            # get positions of real SNPs if any are exist for segment
            if np.any(nslice > 2):
                # get SNP data at position
                isnp = pslice[nslice > 2] + 1
                idx = np.where(np.in1d(snpcnt[0], isnp))[0]
                ref, alt = snpcnt[1][idx], snpcnt[2][idx]
                snp = np.sum(np.array(fnc.count_all_pairs(ssize, alt)), axis=1)
                snp[0] = ns[0] * 0.5 * ssize * (ssize - 1) - snp[1]
            else:
                snp = np.array((ns[0] * 0.5 * ssize * (ssize - 1), 0))
            sdata[ii, 3:] = snp
            # add sim/snp data to lists
            assert np.sum(sim) == np.sum(snp)

        # record exp/obs het around each cons segment in region
        segdata.append(sdata)

    # concatenate all data into single array
    if len(segdata):
        segdata = np.concatenate(segdata)
    else:
        msg = '{} contains no regions below {}'.format(cst.chrom, bth)
        print(msg)

    return segdata


def observed_vs_expected_tail(hom, het, bsx, pred, cnt, bthresh):
    # filter values below bthresh
    bi = (bsx <= bthresh)
    hom, het, bsx, pred, cnt = hom[bi], het[bi], bsx[bi], pred[bi], cnt[bi]

    # get sorting indices based on prediction value
    si = np.argsort(bsx)

    # sort all data segments used for calculations
    hom, het, bsx, pred, cnt = hom[si], het[si], bsx[si], pred[si], cnt[si]

    # simulate data from predictions
    shom, shet = simulated_data(pred, cnt, 216).T

    # calculate pi and simulated pi for each of the bins
    pi = 1.0 * het / (het + hom)
    spi = 1.0 * shet / (shet + shom)

    return np.column_stack((bsx, pi, spi))


def poly_and_conserved_per_morgan(cst):
    """analysis of polymorphisms around conserved regions"""
    assert isinstance(cst, ChromStruct)
    cnts = []
    for ch in cst.chroms:
        # reset chrom
        cst.chrom = ch

        # load genetic map as GeneticMap
        gmp = GeneticMap(ch, cst.gmap_files)

        # create a grid of 0.1cM windows across chromosome
        gstart, gend = gmp.gmap_limits
        wins = np.arange(gstart, gend, 0.1)

        # count neutral sites per 0.1 cM window
        fz = np.load(cst.neut_masks.replace('pr94.nff2.rmneut', 'pr95.nff_2'))
        n_gpos = gmp.interp_gpos(np.where(fz['neutmask'] > 2)[0])
        n_cnt = np.histogram(n_gpos, wins)[0]

        # count conserved sites per 0.1cM window
        ano = AnnoSegments(gmp, afile=cst.bs_target(cst.bs_annos[0]))
        c_cnt = np.histogram(ano.rsites, wins)[0]

        # save array of neut/cons columns for current chrom
        cnts.append(np.column_stack((n_cnt, c_cnt)))

    return np.array(cnts)


def save_cllh_tail(folder_name, bthresh=0.43):
    # find composite file in folder
    fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1
    f_init = fdir + flist[0]
    # initialize ChromStruct with composite init file
    cst = ChromStruct(chrom='chr1', init=f_init)
    cst.vars.num_cores = 1
    ipt = prepare_inputs(cst)
    pred, bsx = get_pred(ipt), get_bsx(ipt)
    cnt = []
    for ch in cst.chroms:
        cst.chrom = ch
        cnt.append(np.load(cst.dv_files)['dv'][:, 0])
    cnt = np.concatenate(cnt)
    cnt = cnt[cnt > 0]
    garr = observed_vs_expected_tail(ipt.hom, ipt.het, bsx, pred, cnt, bthresh)
    fsave = '/'.join(f_init.split('/')[:-1]) + '/het_tail_{}.txt'.format(bthresh)
    np.savetxt(fsave, garr)


def save_cllh_gap(folder_name, num=100):
    # find composite file in folder
    fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1
    f_init = fdir + flist[0]
    # initialize ChromStruct with composite init file
    cst = ChromStruct(chrom='chr1', init=f_init)
    cst.vars.num_cores = 1
    ipt = prepare_inputs(cst)
    pred, bsx = get_pred(ipt), get_bsx(ipt)
    cnt = []
    for ch in cst.chroms:
        cst.chrom = ch
        cnt.append(np.load(cst.dv_files)['dv'][:, 0])
    cnt = np.concatenate(cnt)
    cnt = cnt[cnt > 0]
    garr = observed_vs_expected_cllh(ipt.hom, ipt.het, bsx, pred, cnt, num)
    fsave = '/'.join(f_init.split('/')[:-1]) + '/cllh_gap_{}.txt'.format(num)
    np.savetxt(fsave, garr)


def save_obex_cons(folder_name, bthresh=0.65, param_folder=None):
    f_init = get_comp_file(folder_name)
    # initialize ChromStruct with composite init file
    cst = ChromStruct(chrom='chr1', init=f_init)
    # if params from outside folder used, load
    if param_folder is not None:
        xf_init = get_comp_file(param_folder)
        xparams = ChromStruct('chr1', init=xf_init).params
    else:
        xparams = None
    # get obs/exp data from each segment on each chrom
    seg_data = []
    for ch in cst.chroms:
        cst.chrom = ch
        sdta = obsd_vs_expt_tail_cons(cst, bthresh, xparams)
        if len(sdta):
            seg_data.append(sdta)
    seg_data = np.concatenate(seg_data)
    if xparams is None:
        fsv = '/'.join(f_init.split('/')[:-1]) + \
              '/obex_{}.txt'.format(bthresh)
    else:
        fsv = '/'.join(f_init.split('/')[:-1]) + \
              '/xparam.obex_{}.txt'.format(bthresh)
    np.savetxt(fsv, seg_data)


def main():
    if len(argv) == 2:
        save_obex_cons(argv[1], 0.65)
    elif len(argv) == 3:
        save_obex_cons(argv[1], float(argv[2]))
    elif len(argv) == 4:
        save_obex_cons(argv[1], float(argv[2]), argv[3])
    else:
        print('usage: notebook_2 <folder_name> [b_thresh] [param_folder]')
        exit(1)


if __name__ == '__main__':
    if os.getcwd().startswith('/ifs/'):
        main()
    else:
        print('not called from cluster...exiting')
