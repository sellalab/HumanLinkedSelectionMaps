#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:09:14 2018

@author: davidmurphy
"""

import re
import os
import sys

import numpy as np
from sys import argv
from collections import namedtuple
import data_processing.data_tools as dtl
import data_processing.precalc_tools as prc
from classes.geneticmap import GeneticMap
from classes.runstruct import ChromStruct, root_dir, chromosome_length, \
    human_autosomes


# namedtuple for optimization inputs
OptimizerInputs = namedtuple('OptimizerInputs', 'cst hom het u bs cs s')

# TODO: this should include CS files and perhaps sampling mask as well?
# named tuple for compressed array file paths
ArrayFiles = namedtuple('ArrayFiles', 'f_sg f_bs f_nu f_nt f_dv f_pl')


def build_neutmask(cst):
    """create a neutral mask for the ChromStruct dataset"""

    # phastcons 100-vert p=0 sites
    consmask = dtl.conservation_mask(cst.wigz_ncons, cst.chrom)

    # strict calling mask from 1000 genomes project
    callmask = dtl.sequence_mask(cst.call_mask, cst.chrom)

    # aligned positions to outgroup mask
    submask = dtl.substitution_mask(cst.axt_outgroup, cst.chrom)

    # aggregate filter regions into a binary mask array
    filtmask = dtl.get_filter_mask(cst.chrom, cst.files.fnff)

    # SNP sites (used to encode polymorphism in the final neutrality mask
    print('READING SNP DATA FROM FILE: {}'.format(cst.snp_files))

    # remove CpG SNPs and CpG islands
    if 'filter.CpG' in cst.tkn:
        print('FILTERING CpG SITES AND CpG ISLANDS')
        # # get ancestral state of site to call CpG mutations
        # anc_fname = cst.ancs_files
        # cpg = dtl.cpg_mask(cst.chrom, anc_fname)
        # filtmask[cpg] = False
        # filtmask[cpg + 1] = False
        # filter polymorphic sites from CpGs
        snpcnt = dtl.snpcount(cst.snp_files, returnbases=True)
        issnp, ref, alt = snpcnt[0], snpcnt[1], snpcnt[3]
        pairs = np.core.defchararray.add(ref, alt)
        icpg = np.in1d(pairs, ['GA', 'CT'])
        filtmask[(issnp - 1)[icpg]] = 0
        # add CpG islands to filter mask
        cpgi = np.loadtxt(cst.cpg_isld, usecols=(1, 2), dtype=int)
        filtmask = dtl.mask_segments(filtmask, cpgi)
    else:
        issnp = dtl.snpcount(cst.snp_files)[0]

    # # remove chrX PAR (source: http://genome.ucsc.edu/cgi-bin/hgGateway)
    # if cst.chrom == 'chrX':
    #     filtmask[60000:2699520] = False
    #     filtmask[154931043:155260560] = False

    # neutral mask args (include out file)
    nmsk_args = (consmask, callmask, submask, issnp, filtmask, cst.neut_masks)
    # save and return mask
    nmsk = dtl.neutrality_mask(*nmsk_args)

    return nmsk


def build_neutmask_2(cst):
    """get neutral mask based on conservation, sequence calling and filters"""
    # based on the non conserved choice, find lowest 35% cutoff pmax score
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cst.ncon, c=cst.chrom)
    fdist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=cst.ncon)
    pmax = dtl.dist_percentile(fdist, cst.npct)
    # print 'CALLING NEUTRAL FOR BOTTOM {}% OF {} FILE.'.format(100*cst.npct, fwigz)
    # print 'MAXIMUM PHASTCONS CUTOFF SCORE: {}'.format(pmax)

    # phastcons 100-vert p=0 sites
    if 'cutoff_nval' in cst.tkn:
        wigz_array = dtl.wigz_array(cst.chrom, fwigz, use_neg1=True)
        ncmask = (0 <= wigz_array) & (wigz_array <= cst.nval)
        print('TAKING NEUTRAL SITES BELOW PHASTCONS CUTOFF SCORE: {}'.format(cst.nval))

    else:
        print('CALLING NEUTRAL FOR BOTTOM {}% OF {} FILE.'.format(
            100 * cst.npct, fwigz))
        print('MAXIMUM PHASTCONS CUTOFF SCORE: {}'.format(pmax))
        ncmask = dtl.conservation_mask(fwigz, cst.chrom, pmax=pmax)
    msg = 'NEUTRAL PASSING SCORE COUNT: {:.3e}. FRACTION OF CHROM: {:.1f}%'
    nsum = ncmask.sum()
    print(msg.format(nsum, 100.0 * nsum / cst.chlen))

    # strict calling mask from 1000 genomes project
    callmask = dtl.sequence_mask(cst.call_mask, cst.chrom)
    # aggregate filter regions into a binary mask array
    filtmask = dtl.get_filter_mask(cst.chrom, cst.files.fnff)
    # load genic merge mask separately
    fgene = cst.data + '/coords/genicMerge/{}.genicMerge.bed'.format(cst.chrom)
    gene_segs = np.loadtxt(fgene, dtype=int, usecols=(1, 2))
    gene_mask = dtl.mask_segments(np.ones(shape=ncmask.size, dtype=bool),
                                  gene_segs)

    # remove CpG SNPs and CpG islands
    if 'filter.CpG' in cst.tkn:
        print('FILTERING CpG SITES AND CpG ISLANDS')
        # # get C positions of all CpGs from ancestor
        # icpg = dtl.cpg_mask(cst.chrom, cst.ancs_files)
        # get C positions of all CpGs from hg19 reference
        icpg = dtl.cpg_mask(cst.chrom, cst.refg_files)
        # get snp pos, ref, alt
        snpcnt = dtl.snpcount(cst.snp_files, returnbases=True)
        issnp, ref, alt = snpcnt[0], snpcnt[1], snpcnt[3]
        # create SNP pairs of ref/alt
        pairs = np.core.defchararray.add(ref, alt)
        # find CpG C positions overlapping C>T SNPs
        cti = np.in1d(icpg + 1, issnp[pairs == 'CT'])
        # find CpG G positions overlapping G>A SNPs
        gai = np.in1d(icpg + 2, issnp[pairs == 'GA'])

        # filter each category of CpG SNP
        filtmask[icpg[(cti | gai)]] = False

        # add CpG islands to filter mask
        cpgi = np.loadtxt(cst.cpg_isld, usecols=(1, 2), dtype=int)
        filtmask = dtl.mask_segments(filtmask, cpgi)

    # remove BGC type mutations
    if 'filter.BGC' in cst.tkn:
        is_bgc = ['AC', 'AG', 'CA', 'CT', 'GA', 'GT', 'TC', 'TG']
        print('FILTERING BIASED GENE CONVERSION SNPS')
        # get snp pos, ref, alt
        snpcnt = dtl.snpcount(cst.snp_files, returnbases=True)
        issnp, ref, alt = snpcnt[0], snpcnt[1], snpcnt[3]
        # adjust SNP position to 0 based for masking
        issnp -= 1
        # create SNP pairs of ref/alt
        pairs = np.core.defchararray.add(ref, alt)
        # get a mask of BGC pairs
        ibgc = np.in1d(pairs, is_bgc)
        # mask out BGC pairs in the filter mask
        filtmask[issnp[ibgc]] = False

    # remove CG hypermutable regions
    if 'filter.CG.mut' in cst.tkn:
        print('FILTERING CG HYPERMUTABLE REGIONS')
        # get dict of hypermutable segments for each chrom
        cg_dict = dtl.cg_hypermutable()
        # if current chrom includes hypermutable segments, filter them
        if cst.chrom in cg_dict:
            segs = cg_dict[cst.chrom]
            for (si, sj) in segs:
                filtmask[si:sj] = False
        # print 'INTERSECTING WITH 8 ALIGNABLE SPECIES MASK'
        # f_name = '{}.euarchontoglires.0.35.aligned.8.nmsk.npz'
        # fm8 = root_dir + '/data/nmsk/' + f_name.format(cst.chrom)
        # a8_mask = np.load(fm8)['neutmask']
        # filtmask &= a8_mask

    # initialize each position to True (i.e. neutral)
    neutmask = np.ones(shape=ncmask.size, dtype=bool)
    # remove non-passing sites in polymorphism call mask (bool)
    neutmask &= callmask
    # apply optional filtering mask
    neutmask &= filtmask
    filter_sum = np.sum(neutmask)
    # remove genic merge segments
    neutmask &= gene_mask
    gene_sum = np.sum(neutmask)
    # default mask: 100-vertebrates phastCons score == 0
    neutmask &= ncmask
    neut_sum = np.sum(neutmask)

    # # remove anything that might be misspecified by being on the edge of gmap
    # if 'filter.gmap.edge' in cst.tkn:
    #     print 'FILTERING GMAP EDGES + {}cM'.format(cst.gmsk)
    #     gmsk = np.ones(shape=cst.chlen, dtype=bool)
    #     # load the gmap
    #     gmp = GeneticMap(cst.chrom, cst.gmap_files)
    #     # get genetic map position to mask up to from 5' direction
    #     gpos5 = gmp.interp_gpos(gmp.viable_5prime) + cst.gmsk
    #     # get physical position corresponding to the genetic map masking
    #     pos5 = gmp.interp_pos(gpos5)
    #     # get the index to maximum non-zero position, set region to False
    #     gmsk[:pos5] = False
    #     # get genetic map position to mask starting in the 3' direction
    #     gpos3 = gmp.interp_gpos(gmp.viable_3prime) - cst.gmsk
    #     # get physical position corresponding to the genetic map masking
    #     pos3 = gmp.interp_pos(gpos3)
    #     # get the index to maximum non-zero position, set region to False
    #     gmsk[pos3:] = False
    #     neutmask &= gmsk
    #     gmsk_sum = np.sum(neutmask)
    #     # print sums after each filter step
    #     msg = 'filter {}\ngenic {}\nneutral {}\ngmask {}\n'
    #     print msg.format(filter_sum, gene_sum, neut_sum, gmsk_sum)
    # else:
    #     # print sums after each filter step
    #     msg = 'filter {}\ngenic {}\nneutral {}\n'
    #     print msg.format(filter_sum, gene_sum, neut_sum)

    # remove genetic map edges
    print('FILTERING GMAP EDGES + {}cM'.format(cst.gmsk))
    gmsk = np.ones(shape=cst.chlen, dtype=bool)
    # load the gmap
    gmp = GeneticMap(cst.chrom, cst.gmap_files)
    # get genetic map position to mask up to from 5' direction
    gpos5 = gmp.interp_gpos(gmp.viable_5prime) + cst.gmsk
    # get physical position corresponding to the genetic map masking
    pos5 = gmp.interp_pos(gpos5)
    # get the index to maximum non-zero position, set region to False
    gmsk[:pos5] = False
    # get genetic map position to mask starting in the 3' direction
    gpos3 = gmp.interp_gpos(gmp.viable_3prime) - cst.gmsk
    # get physical position corresponding to the genetic map masking
    pos3 = gmp.interp_pos(gpos3)
    # get the index to maximum non-zero position, set region to False
    gmsk[pos3:] = False
    neutmask &= gmsk
    gmsk_sum = np.sum(neutmask)
    # print sums after each filter step
    msg = 'filter {}\ngenic {}\nneutral {}\ngmask {}\n'
    print(msg.format(filter_sum, gene_sum, neut_sum, gmsk_sum))

    # save to filename if mask file does not exist
    if not os.path.isfile(cst.neut_masks):
        np.savez_compressed(cst.neut_masks, neutmask=neutmask)

    return neutmask


def get_arrays(cst, plen=None):
    """
    Get the basic arrays used for CLLH optimization and results analysis.

    Use reference data and precalculated LS maps. The only synthetic DATA
    file is the neutral mask file, which is generated using PHASTCONS scores,
    callability masks from neutral variation data,
    :param cst:
    :param b_anno:
    :param c_anno:
    :param plen:
    :return:
    """
    # check wheter neut mask file exists, otherwise build it
    if os.path.isfile(cst.neut_masks):
        print('LOADING EXISTING NEUTRAL MASK: {}'.format(cst.neut_masks))
        nmsk = np.load(cst.neut_masks)['neutmask']
    else:
        # TODO: Note to Allen: we will need to prepare all of the intermediate files to rebuild these
        print('NO EXISTING NEUTRAL MASK FOUND. PROVIDE PRECURSOR FILES AND THEN UNCOMMENT CODE BELOW')
        exit(1)
        # build and save, then return
        # print('BUILDING NEUTRAL MASK: {}'.format(cst.neut_masks))
        # nmsk = build_neutmask_2(cst)

    # # create neutmask by default (even if one exists)
    # if 'aligned.8' in cst.tkn:
    #     print('LOADING ALIGNED 8 NEUTRAL MASK: {}'.format(cst.neut_masks))
    #     nmsk = np.load(cst.neut_masks)['neutmask']
    # else:
    #     print('BUILDING NEUTRAL MASK: {}'.format(cst.neut_masks))
    #     nmsk = build_neutmask_2(cst)

    # generate input arrays from reference maps and polymorphism data
    sgarr, bsarr, csarr, nuarr = join_maps(cst, plen, nmsk)
    ntarr, dvarr, plarr = join_poly(cst, nmsk, sgarr)

    return sgarr, bsarr, csarr, nuarr, ntarr, dvarr, plarr


def join_poly(cst, nmsk, sarr):
    """
    Return [hom, het], [site, subs] and [poly_count] for each segment.

    :type cst: ChromStruct
    :param nmsk: neutral mask archive
    :type nmsk: np.lib.npyio.NpzFile
    :param sarr: segments array
    :type sarr: np.ndarray
    :return: nt_arr, dv_arr, ply_arr
    :rtype: np.ndarray
    """
    # get array of snp_pos, ref_count, alt_count from snp files
    print('READING SNP DATA FROM FILE: {}'.format(cst.snp_files))
    snps = np.column_stack(dtl.snpcount(cst.snp_files))

    # get polymorphism & divergence summaries compressed into segments
    pa, sample = prc.compress_data(nmsk, snps, sarr)
    print('SAMPLE SIZE INFERRED FROM SNP FILE: {}'.format(sample))

    # retrieve sample size from SNP data and save to cst if not already set
    # only save for chr1 to avoid colliding open/write/save operations
    if (cst.stat.indv != sample) and cst.chrom == 'chr1':
        cst.stat.indv = sample
        cst.save()

    # return polymorphism, divergence data arranged on segment grid
    return pa[:2].T, pa[2:4].T, pa[4]
    

def load_maps(cst, pl, nmsk):
    """
    Return combined segments and values from precalc maps and neutral u map.

    :param cst:
    :param b_anno:
    :param pl:
    :param nmsk:
    :return:
    """
    # load bmap segment lengths and bvals from each map into separate lists

    # collect bmap files
    mrg = True if pl else False
    bmap_files = []
    bmsg = 'READING SEGMENTS AND VALUES FROM BMAP: {}'
    for (ba, tvec) in zip(cst.bs_annos, cst.bdfe):
        for ti in tvec:
            f_bm = cst.bkgd_file(ba, ti, plen=pl, merge=mrg)
            print(bmsg.format(f_bm))
            bmap_files.append(f_bm)
    # collect bmap segments and values
    bsegs, bvals = prc.collect_lsmaps(cst.chrom, bmap_files)

    # collect cmap files
    cmap_files = []
    cmsg = 'READING SEGMENTS AND VALUES FROM CMAP: {}'
    for (ca, svec) in zip(cst.cs_annos, cst.cdfe):
        for si in svec:
            f_cm = cst.cmap_file(ca, si)
            print(cmsg.format(f_cm))
            cmap_files.append(f_cm)
    # collect cmap segments and values
    csegs, cvals = prc.collect_lsmaps(cst.chrom, cmap_files)

    # # load segment lengths and rates for neutral sub maps
    # sub_counts = dtl.substitution_counts(nmsk, 2e4, cst.chrom)
    # usegs, uvals = prc.prepare_umap(cst.chrom, sub_counts)

    # # load new segment rates from files
    # fsub = '{}/phast/sub_rates/{}.sub_rates.npz'.format(cst.data, cst.chrom)
    # if 'filter.CpG' in cst.tkn:
    #     print 'USING CpG ADJUSTED SUBSTITUTION RATES'
    #     subrate = np.load(fsub)['rates'][:, (0, 2)]
    # else:
    #     print 'USING TOTAL SUBSTITUTION RATES'
    #     subrate = np.load(fsub)['rates'][:, :2]
    # usegs, uvals = prc.prepare_umap_2(cst.chrom, subrate)

    # load new segment rates from files
    # fsub = '{}/phast/sub_rates/{}.sub_counts.npz'.format(cst.data, cst.chrom)

    # /ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs/
    # aligned_8_win_10000_subcounts/chr1.subcount.filled.txt
    # rdir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/runs/'
    # sdir = rdir + 'aligned_{}_win_{}'.format(cst.nspc, cst.wind)
    # sdir = rdir + 'aligned{}_win{}'.format(cst.nspc, cst.wind)
    # if cst.slid:
    #     pct = int(100.0 * 1.0 / cst.slid)
    #     #     fout = '{}/{}.subcount.filled.{}pct.txt'.format(pth, ch, pct)
    #     fsub = '{}/{}.subcount.filled.{}pct.txt'.format(sdir, cst.chrom, pct)
    #     print('LOADING SLIDING SUBSTITUTION COUNT FILE {}'.format(fsub))
    # else:
    #     fsub = '{}/{}.subcount.filled.txt'.format(sdir, cst.chrom)
    #     print('LOADING NON-OVERLAPPING SUBSTITUTION COUNT FILE {}'.format(fsub))
    # if 'filter.CpG' in cst.tkn:
    #     print('ESTIMATING NU USING CpG ADJUSTED SUBSTITUTION COUNTS')
    #     subcnt = np.loadtxt(fsub)[:, (0, 1, 2, 4)]
    # elif 'filter.BGC' in cst.tkn:
    #     print('ESTIMATING NU USING BGC ADJUSTED SUBSTITUTION COUNTS')
    #     subcnt = np.loadtxt(fsub)[:, (0, 1, 2, 5)]
    # else:
    #     print('ESTIMATING NU USING TOTAL SUBSTITUTION COUNTS')
    #     subcnt = np.loadtxt(fsub)[:, :4]

    print('ESTIMATING NU USING TOTAL SUBSTITUTION COUNTS')
    subcnt = np.loadtxt(cst.neut_subs)[:, :4]
    usegs, uvals = prc.prepare_umap_4(cst.chrom, subcnt)

    # return combined lists of segments and values
    segs = bsegs + csegs + [usegs]
    vals = bvals + cvals + [uvals]
    return segs, vals


def join_maps(cst, pl, nmsk):
    # load all segs and vals
    segs, vals = load_maps(cst, pl, nmsk)

    # combine data into unified segment grid
    sgarr, v = prc.combine_segments(cst.chrom, segs, vals)

    # split segmented data using param indices
    bvals = v[:, cst.fixed.bi:cst.fixed.bj]
    cvals = v[:, cst.fixed.ci:cst.fixed.cj]
    uvals = v[:, cst.fixed.cj:cst.fixed.cj+1]

    # return segment grid, bvals cvals and uvals arranged on segment grid
    return sgarr, bvals, cvals, uvals


def load_saved(cst, chroms=None):
    """
    Load saved LH input arrays.

    :param cst:
    :return: nt, u, bs, d, seg UN-masked
    :rtype: tuple
    """
    # raise flag if segment count stats haven't been added to cst
    if (cst.stat.totsegs == []) or (cst.stat.msksegs == []):
        # blank both lists
        cst.stat.totsegs = []
        cst.stat.msksegs = []
        count_segs = True
    else:
        count_segs = False

    if chroms is not None:
        count_segs = False
    else:
        chroms = cst.chroms

    sg, bs, cs, nu, nt, dv, pl = [], [], [], [], [], [], []
    for ch in chroms:
        # open npz archive of all input arrays
        cst.chrom = ch

        # if skip chrom is set to the current chrom, don't load it
        if ch == cst.vars.drop_chrom:
            print('DROPPING CHROMOSOME {} FROM DATASET'.format(ch))
            continue

        # count neutral sites per segment
        div = np.load(cst.dv_files)['dv']
        ntsites = div[:, 0]

        # count segments and segments with DATA
        if count_segs:
            cst.stat.totsegs.append(ntsites.size)
            cst.stat.msksegs.append(np.sum(ntsites > 0))

        # collect arrays into respective lists
        dv.append(div)
        sg.append(np.load(cst.sg_files)['sg'])
        nu.append(np.load(cst.nu_files)['nu'])
        nt.append(np.load(cst.nt_files)['nt'])
        pl.append(np.load(cst.pl_files)['pl'])

        # check if bs is present in the run
        if cst.bnum > 0:
            bs.append(np.load(cst.bs_files)['bs'])
        else:
            bs.append([])

        # check if cs is present in the run
        if cst.cnum > 0:
            cs.append(np.load(cst.cs_files)['cs'])
        else:
            cs.append([])

    # # if segment counts stats were collected, save the cst file
    # if count_segs:
    #     cst.save()

    # concatenate arrays across chromosomes
    sg = np.concatenate(sg)
    bs = np.concatenate(bs)
    cs = np.concatenate(cs)
    nu = np.concatenate(nu)
    nt = np.concatenate(nt)
    dv = np.concatenate(dv)
    pl = np.concatenate(pl)

    # set empty LS maps to None
    if not bs.size:
        bs = None
    if not cs.size:
        cs = None

    # generate a jackknife mask if jackknife is in use
    if cst.vars.use_jackknife:
        jidx, jwin = cst.vars.jackknife_index, cst.vars.jackknife_window
        print('USING JACKKNIFE SAMPLING: INDEX={} WINDOW={}'.format(jidx, jwin))
        m = get_jackknife_mask(sg, jidx, jwin)
        # if jackknife mask has no effect on data, stop program
        n_1, n_2 = dv.sum(), dv[m].sum()
        if n_1 == n_2:
            msg = 'ERROR: jackknife index {} and window {} have no effect.'
            print(msg.format(jidx, jwin))
            exit(1)
        else:
            d = n_1 - n_2
            print('JACKKNIFE MASKING REMOVED {} OF {} SITES'.format(d, n_1))
        # otherwise just return the masked data
        sg, nu, nt, dv, pl = sg[m], nu[m], nt[m], dv[m], pl[m]
        # check that bs/cs arrays exist before masking
        if bs is not None:
            bs = bs[m]
        if cs is not None:
            cs = cs[m]

    return sg, bs, cs, nu, nt, dv, pl


def adjust_arrays(cst, bs, cs, nu, nt, dv, pl):
    """mask empty segments and apply bounds precalc maps if specified"""
    # take site counts from dv, mask filters segments with no neutral sites
    dv = dv[:, 0]
    msk = dv > 0

    # TODO: why is this step suddenly needed if array hasn't changed shape?
    if len(nu.shape) == 2:
        nu = nu[:, 0]

    # get a downsampling mask if specified
    if cst.vars.down_sample:
        ds_msk = get_dsmask(cst, msk.sum())
    else:
        ds_msk = None

    # mask and aggregate all the arrays
    adjusted = []
    for arr in (bs, cs, nu, nt, dv, pl):
        # mask out blocks without neutral data
        if (arr is not None) and (ds_msk is None):
            adjusted.append(arr[msk])
        # mask out blocks without neutral data and downsample
        elif (arr is not None) and (ds_msk is not None):
            adjusted.append(arr[msk][ds_msk])
        # record "None" for empty array
        else:
            adjusted.append(None)

    bs, cs, nu, nt, dv, pl = adjusted

    # convert nu to a constant if set in variables
    if cst.vars.mu_const:
        assert cst.stat.meandiv
        sys.stderr.write('USING CONSTANT MU = {}\n'.format(cst.stat.meandiv))
        sys.stdout.flush()
        nu = np.full(shape=len(dv), fill_value=cst.stat.meandiv)

    # determine whether map is old or new by inspection
    if bs is not None:
        # USE BY DEFAULT!
        bs *= np.log1p(-1.0 / cst.bscl)
        #
        # # only observed in new maps
        # if bs.max() > cst.bscl:
        #     # rescale binned little b values
        #     bs *= np.log1p(-1.0 / cst.bscl)
        # else:
        #     # convert binned big B to little b with min log(0.01)
        #     assert np.unique(bs).size <= cst.bscl + 1

        # bound bs by preset minimum if specified and bs is present
        if cst.fixed.min_bs is not None:
            bs = np.maximum(cst.fixed.min_bs, bs)

    return bs, cs, nu, nt, dv, pl, msk


def get_jackknife_mask(segs, jidx, jwin):
    """get mask to remove jackknife segment from data"""
    # get the start/end POSITION based on index and size
    i, j = jidx*jwin, (jidx+1)*jwin
    # convert segment lengths to positions
    pos = np.cumsum(segs)
    # check that segment sum makes sense
    assert pos[-1] == sum(chromosome_length(c) for c in human_autosomes)
    # find the start/end INDICES in the data to remove jackknife slice
    ji, jj = np.searchsorted(pos, [i, j])
    # check that ji makes sense
    assert ji <= len(pos) - 1
    # bound the maximum jj value
    jj = min(len(pos), jj)
    # create all TRUE mask, flip jackknife positions to FALSE
    jmsk = np.ones(shape=len(pos), dtype=bool)
    jmsk[ji:jj] = False

    return jmsk


def get_dsmask(cst, size, makenew=False):
    """creates a new downsampling mask or loads existing"""
    # load existing mask
    if os.path.isfile(cst.dsmsk_file) and not makenew:
        mask = np.load(cst.dsmsk_file)['dsmsk']
    # save new mask to file name and return
    else:
        p = cst.vars.down_sample_pct
        mask = downsampling_mask(size=size, p=p)
        np.savez_compressed(cst.dsmsk_file, dsmsk=mask)

    return mask


# TODO: should sampling mask "s" be carried through here?
def prepare_inputs(cst):
    """
    Concatenate arrays and group; split into chunks if num_cores > 1.

    :param cst: ChromStruct containing dirs, params and constants
    :param new: flag indicating whether file dir structure is old or new
    :param loaded_arrays: tuple of arrays returned by "load_saved" function
                          (optional; default=None)
    :return: set of inputs for optimize program
    :type cst: ChromStruct
    :type new: bool
    :type loaded_arrays: tuple
    :rtype: container
    """
    # load compressed arrays
    saved_arrs = load_saved(cst)[1:]

    # mask empty segments, modify arrays if specified (don't return mask)
    bs, cs, nu, nt, dv = adjust_arrays(cst, *saved_arrs)[:5]

    # reshape nt array so hom & het are memory neutral views on nt array
    nt = nt.T
    hom, het = nt

    # take basic pre-inference stats
    # print 'ARRAY SHAPES: dv={}, nu={}'.format(dv.shape, nu.shape)
    cst = summary_stats(cst, hom, het, nu, dv)

    # organize inputs into single arrays and return
    inputs = OptimizerInputs(cst, hom, het, nu, bs, cs, None)

    # split inputs for multi-threading if more than one core in use
    if cst.vars.num_cores > 1:
        inputs = split_input(inputs)

    return inputs


def split_input(inputs):
    # get indices to break up data into evenly sized chunks
    ij = split_index(inputs.cst.stat.used_segs, inputs.cst.vars.num_cores)

    # divide input arrays into slices using the split indices
    split_inputs = []
    for (i, j) in ij:
        # every input slice starts with the ChromStruct
        input_slice = [inputs.cst]
        for a in inputs[1:]:
            # add "None" to unfilled input fields
            if a is None:
                input_slice.append(None)
            # add sliced array for filled input fields
            else:
                input_slice.append(a[i:j])

        # add list of slices to list of chunks
        split_inputs.append(OptimizerInputs._make(input_slice))

    # clear old inputs from memory
    del inputs

    # return list of input chunks
    return split_inputs


def split_index(n_sites, m_blocks):
    """
    Get indices to split n_sites into m_blocks of data.

    A function that returns slice indices (i and j) that divide data of n_sites
    into m_blocks of equal size, except for the final chunk, which spans the
    remaining data and may include up to m_blocks-1 additional data points.
    NOTE: when n_sites >> m_blocks then step >> m_blocks and the added data in
    the final chunk is negligible.

    :param n_sites: number of segments in each data array
    :param m_blocks: number of chunks data arrays should be broken into
    :return i, j: (lower, upper) slice indices, respectively
    """
    # step size that indexes data into n chunks
    step = int(n_sites / m_blocks)
    # remaining sites if n_sites not evenly divisible by m_blocks
    remain = n_sites % m_blocks
    # upper edge chunk indices
    j = np.arange(step, n_sites+1, step, dtype='uint')
    # lower edge chunk indices
    i = j - step
    # correct upper edge of the last chunk for remainder
    j[-1] += remain
    # stack [lower, upper] indices
    indices = np.column_stack((i, j))
    # final segment is at MOST m_blocks longer than other segments
    assert (indices[-1, 1] - indices[-1, 0]) < (step + m_blocks)
    # last upper edge index = n_sites
    assert j[-1] == n_sites

    return indices


def summary_stats(cst, hom, het, u, d):
    # get essential stats to start optimzation
    cst.stat.used_segs = len(hom)
    cst.stat.hets = het.sum()
    cst.stat.homs = hom.sum()
    cst.stat.pairs = cst.stat.hets + cst.stat.homs
    cst.stat.meanpi = cst.stat.hets / cst.stat.pairs
    cst.stat.meandiv = np.average(u, weights=d, axis=0)
    cst.fixed.tau_init = cst.stat.meandiv / cst.stat.meanpi

    # only set tau if it is in the default initial state
    if cst.params[-1] < 0:
        cst.params[-1] = cst.fixed.tau_init

    return cst


def downsampling_mask(size, p):
    """
    Generate a mask to downsample dataset.

    :param size: size of array to downsample
    :param p: fraction of the array to randomly choose
    :return: boolean mask for blocked dataset
    :type size: int
    :type p: float
    :rtype: np.ndarray
    """
    # expected fraction of "False" cells in the mask
    q = 1.0 - p
    # generate mask using weighted random choice of True/False
    mask = np.random.choice([True, False], size=size, p=[p, q])

    return mask


def save_arrays(cst, plen):
    """
    Get compressed data and precalc map arrays.

    :type cst: ChromStruct
    :param b_anno: annotation for BS files
    :param plen:
    :rtype: None
    """
    # make directories for saving compressed data files:
    cmpdir = root_dir + '/compress'
    os.makedirs(cmpdir, exist_ok=True)
    for d in ['parr', 'darr', 'narr', 'uarr', 'barr', 'carr', 'sarr']:
        os.makedirs('{}/{}'.format(cmpdir, d), exist_ok=True)

    # get arrays of data compressed into segments:
    # sg: segment arrays
    # bs: background selection map arrays
    # cs: classic sweep map arrays
    # nu: neutral mutation rate arrays
    # nt: neutral polymorphism data arrays
    # dv: divergence data arrays (not used in current build)
    # pl: polymorphic site count arrays (not used in current build)
    sg, bs, cs, nu, nt, dv, pl = get_arrays(cst, plen)

    # save each compressed data array
    np.savez_compressed(cst.sg_files, sg=sg)
    if len(cst.bs_annos):
        np.savez_compressed(cst.bs_files, bs=bs)
    if len(cst.cs_annos):
        np.savez_compressed(cst.cs_files, cs=cs)
    np.savez_compressed(cst.nu_files, nu=nu)
    np.savez_compressed(cst.nt_files, nt=nt)
    np.savez_compressed(cst.dv_files, dv=dv)
    np.savez_compressed(cst.pl_files, pl=pl)

    return None


def main():
    # # local
    # if root_dir.startswith('/Users/MURPHYD/'):
    init = root_dir + '/result/init_files/YRI.cadd94_gmask_v1.6_without_bstat.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    # this is the size of the chunks used to make B maps in fragments
    plen = 1e7
    # chrom = 'chr22'
    for chrom in human_autosomes:
        # initialize ChromStruct from file
        cst = ChromStruct(chrom=chrom, init=init)
        # generate all the compressed arrays and save
        print('Preparing inputs for {}...\n'.format(chrom))
        save_arrays(cst, plen)
        print('\nDone preparing inputs for {}.'.format(chrom))

    # # remote
    # else:
    #     if len(argv) != 4:
    #         print('usage: lh_inputs <chrom> <init> <plen>')
    #         exit()
    #     chrom = argv[1]
    #     init = argv[2]
    #     plen = eval(argv[3])
    #
    # # initialize ChromStruct from file
    # cst = ChromStruct(chrom=chrom, init=init)
    #
    # # generate all the compressed arrays and save
    # save_arrays(cst, plen)
    # print('\nDone preparing inputs for {}.'.format(chrom))


if __name__ == '__main__':
    main()
