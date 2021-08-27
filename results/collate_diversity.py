
import os
import numpy as np
from sys import argv, stderr, stdout
import data_processing.data_tools as dtl
from data_processing.functions import calc_pi
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import load_saved, adjust_arrays
from classes.runstruct import RunStruct, ChromStruct, root_dir, cst_from_fldr
from classes.geneticmap import GeneticMap
from classes.substitution import recursive_search


__author__ = 'davidmurphy'


def get_yri_params(fldr):
    """get just the predictions for YRI map"""
    cst = cst_from_fldr(fldr)
    return cst.params, cst.stat.meanpi


def err_msg(errstring):
    stderr.write(errstring + '\n')
    stdout.flush()


def main_remote():
    # if len(argv) != 4:
    #     print 'usage: collate_diversity <folder> <width> <focal>'
    #     exit(1)

    if len(argv) == 4:
        # get init_file path and collated bin width from command line args
        # init_file = argv[1]
        width = eval(argv[2])
        focal = argv[3]

        fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
        flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
        assert len(flist) == 1

        # create paths to rsq and final init file
        f_init = fdir + flist[0]

        # initialize RunStruct, adjust collated plot variables if specified
        cst = ChromStruct(chrom='chr1', init=f_init)
        # cst = cst_from_fldr(argv[1])
        if width is not None:
            cst.vars.collated_bin = width
        if focal is not None:
            cst.fcl = focal
        # cst.save()

        # calculated sum of divergence, pi, prediction, sites across bins
        collate_array = collate_diversity(cst)

        # construct out-file
        fpth = '/'.join(f_init.split('/')[:-1])
        fmt = fpth + '/{}.collate.{:.2e}width.npy'
        fout = fmt.format(focal, width)

        # format: bin, pi, div, pred, cnts
        np.save(fout, collate_array)

    elif len(argv) == 5:
        width = eval(argv[2])
        focal = argv[3]
        map_pop = argv[4]
        fdir = root_dir + '/result/final_files/{}/'.format(argv[1])
        flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
        assert len(flist) == 1

        # create paths to rsq and final init file
        f_init = fdir + flist[0]

        # initialize RunStruct, adjust collated plot variables if specified
        cst = ChromStruct(chrom='chr1', init=f_init)
        if width is not None:
            cst.vars.collated_bin = width
        if focal is not None:
            cst.fcl = focal
        # cst.save()

        # calculated sum of divergence, pi, prediction, sites across bins
        collate_array = collate_diversity(cst, map_pop=map_pop)

        # construct out-file
        fpth = '/'.join(f_init.split('/')[:-1])
        fmt = fpth + '/YRI_sorted.{}.collate.{:.2e}width.npy'
        fout = fmt.format(focal, width)

        # format: bin, pi, div, pred, cnts
        np.save(fout, collate_array)
    else:
        print 'usage_1: collate_diversity <folder> <width> <focal>'
        print 'usage_2: collate_diversity <folder> <width> <focal> <map_pop>'
        exit(1)

def main_local():
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS'
    f = rdir + '/result/final_files/YRI.pr95.cleanrun.BS1.6.CS0.0.170714065758.final.txt'
    rst = RunStruct(init=f)
    # rst.vars.collated_bin = 5e-3
    rst.focal = 'syn'

    # gather pi, divergence and prediction summaries across bins
    collate_array = collate_diversity(rst)

    # construct out-file
    out_file = '{root}/result/collate/{label}.{focal}.{neut}.collated.npy'.format(**rst.dict)

    # save collated summaries array to out file
    np.save(file=out_file, arr=collate_array)

    # rst = RunStruct(init=f)
    # in_file = '{rt}/result/collate/{fo}/{lb}.{fo}.{pp}.npy'.format(**rst.dict)
    # # in_file = '{rt}/result/collate/{fo}/{lb}.{fo}.{pp}-1e_3.npy'.format(**rst.dict)
    # x, pi, dv, pr, n = np.load(in_file).T
    # # '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/collate/nsSubs/mcvicker.map.BS1CS0.nsSub.YRI.npy'
    # # '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/collate/nsSub/mcvicker.map.BS1CS0.nsSub.YRI.npy'
    # mcfile = '{rt}/result/collate/nsSub/mcvicker.map.BS1CS0.nsSub.YRI.npy'.format(**rst.dict)
    # mx, mpi, mdv, mpr, mn = np.load(mcfile).T
    #
    # from functions import savitzky_golay as svg
    # import matplotlib.pyplot as plt
    # import seaborn
    #
    # plt.figure(figsize=(11, 7))
    # plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.94, wspace=0.05, hspace=0.1)
    # # plt.subplot(121)
    # window = 35
    # order = 11
    # pi /= (n * rst.stat.meanpi0)
    # pr /= (n * rst.stat.meanpi0)
    # mpi /= (mn * rst.stat.meanpi0)
    # mpr /= (mn * rst.stat.meanpi0)
    # # pi = svg(pi / (n * rst.stat.meanpi0), window_size=window, order=order)
    # # mc = svg(mpr / (mn * rst.stat.meanpi0), window_size=11, order=3)
    # # pr = svg(pr / (n * rst.stat.meanpi0), window_size=21, order=5)
    # plt.plot(x, pi, label='Observed', lw=2.5, color='DarkSlateGray')
    # plt.plot(mx, mpr, label='McVicker\'s prediction', lw=2.5, color='DarkTurquoise', alpha=0.8)
    # plt.plot(x, pr, label='Our prediction', lw=2.5, color='Fuchsia', alpha=0.8)
    #
    # # plt.plot(mx, mpi, label='Observed', lw=2.5, color='MidnightBlue')
    # plt.xlabel('distance to nearest substitution (cM)', fontsize=24)
    # plt.xticks(fontsize=22)
    # plt.xlim(-0.4, 0.4)
    # plt.ylabel('scaled diversity', fontsize=24)
    # plt.yticks(fontsize=22)
    # plt.ylim(0.61, 1)
    # plt.legend(prop={'size': 22}, ncol=1, loc='lower left')
    # # plt.subplot(122)
    # # plt.plot(x, dv / n)
    # plt.show()


def collate_diversity(cst, map_pop=None):
    """
    Collate diversity by genetic distance around focal sites genome wide.
    :param init: RunStruct init file
    """
    # get concatenated positions, predictions and indices per chrom
    # mst = MapStruct(init=init)
    # positions = mst.pos
    # prediction = mst.prediction()
    # chrom_indices = mst.mbnds
    # del mst  # (delete memory-heavy MapStruct)

    # use init file to create DataStruct
    # cst = ChromStruct('chr1', init=init)

    # the bins used to organize predicted and observed diversity data
    span, width = cst.vars.collated_span, cst.vars.collated_bin
    bins = np.arange(-span, span + width, width)

    # store bin coordinates, sums of [pi, pred, counts] per bin for all chroms
    collated_data = np.zeros(shape=[len(bins) - 1, 5])
    collated_data[:, 0] = bins[:-1]

    from datetime import datetime as dt
    total_start = dt.now()
    for ch in cst.chroms:
        start = dt.now()

        # update internal chrom parameter for DataStruct: updates properties
        cst.chrom = ch
        # load complete array data
        sg, bs, cs, nu, nt, dv, pl = load_saved(cst, chroms=[ch])

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

        # convert segments into positions for r r squared calc
        positions = np.cumsum(sg)[msk]

        # get predictions from YRI CADD map
        # TODO: might not want to hardcode this!
        if map_pop is not None:
            err_msg('using YRI CADD parameters for predictions')
            yri_fldr = 'cadd94_gmask_mnb_378'
            # get bmap params from YRI, use theta from specific population
            yri_prm, yri_meanpi = get_yri_params(yri_fldr)
            # get prediction from YRI map
            yri_pred = predicted_pi(yri_prm, cst, nu, bs, cs)
            # rescale predictions to fit the mean pi of current population
            yri_pred *= (cst.stat.meanpi / yri_meanpi)
            # use the YRI CADD predictions as the predictions for collating
            predictions = yri_pred
        else:
            # calculate predicted pi from arrays and average of top 3 params
            predictions = predicted_pi(cst.params, cst, nu, bs, cs)

        # # slice out positions and predictions for the current chromosome
        # bi, bj = chrom_indices[cst.cidx:cst.cidx+2]
        chrom_pred = positions, predictions

        # cM, div, pi, pred for nearest focal within +/- span
        chrom_data = get_distances(cst, chrom_pred)

        # get binned pi, div, pred, cnts per bins defined by span
        collated_data[:, 1:] += bin_distances(cst, chrom_data)

        msg = '{} processing time = {}'.format(ch, dt.now() - start)
        err_msg(msg)

    msg = 'total processing time = {}'.format(dt.now() - total_start)
    err_msg(msg)

    return collated_data


def get_distances(cst, mappred):
    """
    Compress [distance, pairs, hets, pred] data into scale bins
    :param cst: DataStruct for a given chromosome
    :param mappred: [cum_pos, prediction] from the map for the current chrom
    :return pi_obs, pi_pred, count: each bin used (even those without data),
    observed and predicted pi,
    count of sites per bin
    """
    # load neutral mask to get neutral site positions and neutral divergence data (then delete)
    neutmask = np.load(cst.neut_masks)['neutmask']
    neut_pos = np.where(neutmask != 0)[0]

    # create empty neutral array to hold [dist, div, pi, pred] for each neutral site
    narr = np.zeros(shape=(4, len(neut_pos)), dtype='f8')
    narr[1] = neutmask[neut_pos] - 1  # store divergence (adjusted down from 1/2 to 0/1) in column 1
    del neutmask

    # distance to nearest focal, data index for each focal (+1 to neut_pos for 1-based coords to match gmap)
    span = cst.vars.collated_span

    # list the different types of focals that may be used
    point_focals = ['nonsyn', 'YRI_syn_s1', 'hc_derived_cons']
    # segment_focals = ['exon', 'ape_cons94_nonexonic', 'ape_cons94_exonic']
    cons_segments = ['fish_cons94_new_nonexonic', 'fish_cons94_new_exonic',
                     'cadd94_gmask_exonic', 'cadd94_gmask_nonexonic']
    if cst.fcl in point_focals:
        err_msg('using focal points for {}'.format(cst.fcl))
        f_focal = cst.cs_target(cst.fcl)
        if f_focal.endswith('.npy'):
            foc_pos = np.load(f_focal)
        else:
            assert f_focal.endswith('.npz')
            z_focal = np.load(f_focal)
            assert len(z_focal.keys()) == 1
            foc_pos = z_focal[z_focal.keys()[0]]
        dist, idx = distance2focals(neut_pos+1, foc_pos, cst.gmap_files, span)
    elif cst.fcl == 'exon':
        err_msg('using slow function for exons'.format(cst.fcl))
        # assert cst.fcl in segment_focals
        f_focal = cst.bs_target(cst.fcl)
        segments = np.loadtxt(f_focal, usecols=(1,2))
        dist, idx = distance2segments(neut_pos, segments, cst, span)
        # dist, idx = distance2segments_2(neut_pos, segments, cst, span)
    else:
        assert cst.fcl in cons_segments
        err_msg('using conserved segment midpoints for {}'.format(cst.fcl))
        f_focal = cst.bs_target(cst.fcl)
        start, end = np.loadtxt(f_focal, usecols=(1,2)).T
        width = end-start
        mid_pos = start+0.5*width
        dist, idx = distance2focals(neut_pos, mid_pos, cst.gmap_files, span)

    # use argsort to get indices ordering data from from -span to +span along the x-axis
    sidx = dist.argsort()
    dist, idx = dist[sidx], idx[sidx]
    del sidx

    # load neutral poly data for the chrom for the subset indices, calculate pi at each site
    # poly_pos, ref, alt = np.load(dst.neut_poly)['neutpoly'].T
    poly_pos, ref, alt = dtl.snpcount(cst.snp_files)
    # adjust poly_pos to 0 based
    poly_pos -= 1

    # find neutral array indices for poly sites
    nidx = np.where(np.in1d(neut_pos, poly_pos))[0]
    print '{} neutral poly = {}'.format(cst.chrom, nidx.size)

    # there should be fewer poly NEUTRAL than overall poly sites
    assert len(nidx) < len(poly_pos)
    # get the poly index where sites are neutral for ref/alt
    pmsk = np.in1d(poly_pos, neut_pos)
    assert np.sum(pmsk) == len(nidx)

    # set pi values for poly sites at the proper indices in column 2
    narr[2, nidx] = calc_pi(ref[pmsk]+alt[pmsk], alt[pmsk])

    # interpolate prediction for each neutral position
    narr[3] = np.interp(x=neut_pos, xp=mappred[0], fp=mappred[1])  # store pred in column 3

    # filter and sort data using the index found in distance2focals function
    narr = narr[:, idx]

    # finally, set genetic distance in column 0
    narr[0] = dist

    return narr


def distance2focals(neut_pos, foc_pos, gmap_file, span):
    """
    For each neutral position, calculate the genetic distance to the nearest focal position. Return the distances for
    all neutral positions within +/- "span" of a focal position.
    """
    # load genetic map physical and cM columns (these are cols 0 & 2 for our map)
    gmap_bp, gmap_cm = np.loadtxt(gmap_file, usecols=(0, 2)).T

    # convert neutral and focal positions to cM via interpolation on the genetic map
    neut_gpos = np.interp(x=neut_pos, xp=gmap_bp, fp=gmap_cm)
    foc_gpos = np.interp(x=foc_pos, xp=gmap_bp, fp=gmap_cm)

    # indices of nearest focal where focal-cM >= neutral-cM, e.g. nearest focal 3' of neutral sites
    # (use minimum in case neutral site idx is assigned above max focal idx)
    i3prime = np.minimum(np.searchsorted(foc_gpos, neut_gpos), len(foc_gpos)-1)
    # indices of nearest 5' focal, e.g. one step back from the nearest 3' (unless that is 0)
    i5prime = np.maximum(i3prime - 1, 0)

    # get the distances to 5'/3' flanking focal for each neutral site
    gdist_3prime = foc_gpos[i3prime] - neut_gpos
    gdist_5prime = foc_gpos[i5prime] - neut_gpos

    # inequalities should only be violated with neut gpos outside focal gpos range
    assert np.all(gdist_3prime >= 0) or (neut_gpos.max() > foc_gpos.max())
    assert np.all(gdist_5prime <= 0) or (neut_gpos.min() < foc_gpos.min())

    # get the data indices of all positions with distances within the defined range
    idx_3prime = np.where((abs(gdist_3prime) <= abs(gdist_5prime)) & (abs(gdist_3prime) <= span))[0]
    idx_5prime = np.where((abs(gdist_5prime) < abs(gdist_3prime)) & (abs(gdist_5prime) <= span))[0]

    # return the distances and the sorting data indices
    distance = np.concatenate((gdist_5prime[idx_5prime], gdist_3prime[idx_3prime]))
    indices = np.concatenate((idx_5prime, idx_3prime))

    return distance, indices


def distance2segments(neut_pos, segments, cst, span):
    """
    For each neutral position, calculate the genetic distance to the nearest focal position. Return the distances for
    all neutral positions within +/- "span" of a focal position.
    """
    # load genetic map physical and cM columns (these are cols 0 & 2 for our map)
    assert isinstance(cst, ChromStruct)
    gmp = GeneticMap(cst.chrom, cst.gmap_files)

    # sites 5' of segments
    gdist5, gdist3 = [], []
    index5, index3 = [], []

    # special case: first segment
    # get the 5' end of the first segment
    first_seg5 = segments[0,0]
    # get all the neutral sites below this segment
    first_neut = neut_pos[neut_pos < first_seg5]
    # get the genetic map position of the 5' end of the first segment
    first_seg5_gpos = gmp.interp_gpos(first_seg5)
    # get the genetic map position of the neutral sites below this segment
    first_neut_gpos = gmp.interp_gpos(first_neut)
    # get the distances from neutral sites to the first 5' end of first segment
    first_gdist = first_neut_gpos - first_seg5_gpos
    assert np.all(first_gdist <= 0)
    # record the distances to the 5' distances list
    gdist5.append(first_gdist)
    # record the indices to the 5' indices list
    index5.append(np.where(neut_pos < first_seg5)[0])

    # go through one segment at a time and find the neutral sites in between
    # then determine if the 3' end of seg1 or the 5' end of seg2 is closer
    for i in xrange(len(segments)-1):
        # get the 3' end of the first segment and its genetic map position
        pos3 = segments[i,1]
        gpos3 = gmp.interp_gpos(pos3)
        # get the 5' end of the second segment and its genetic map position
        pos5 = segments[i+1,0]
        gpos5 = gmp.interp_gpos(pos5)
        # get the set of indices for neutral sites between the segments
        cur_idx = np.where((neut_pos>pos3) & (neut_pos<pos5))[0]
        # get the subset of neutral sites between the segments and their gpos
        cur_neut = neut_pos[cur_idx]
        cur_gpos = gmp.interp_gpos(cur_neut)
        # get the genetic distances from all neutral sites to 3' and 5' points
        g3 = cur_gpos - gpos3
        assert np.all(g3 >= 0)
        g5 = cur_gpos - gpos5
        assert np.all(g5 <= 0)
        # get mask of sites where g3 is smaller than abs(g5)
        msk = (g3<abs(g5))
        # record the distances and indices to respective lists
        gdist3.append(g3[msk])
        index3.append(cur_idx[msk])
        gdist5.append(g5[~msk])
        index5.append(cur_idx[~msk])

    # special case: last segment
    # get the 3' end of the last segment
    last_seg3 = segments[-1, 1]
    # get all the neutral sites above this segment
    last_neut = neut_pos[neut_pos > last_seg3]
    # get the genetic map position of the 3' end of the last segment
    last_seg3_gpos = gmp.interp_gpos(last_seg3)
    # get the genetic map position of the neutral sites above this segment
    last_neut_gpos = gmp.interp_gpos(last_neut)
    # get the distances from neutral sites to the last 3' end of last segment
    last_gdist = last_neut_gpos - last_seg3_gpos
    # record the distances to the 5' distances list
    gdist3.append(last_gdist)
    # record the indices to the 5' indices list
    index3.append(np.where(neut_pos > last_seg3)[0])

    distance = np.concatenate(gdist3+gdist5)
    indices = np.concatenate(index3+index5)
    span_mask = (abs(distance) <= span)

    return distance[span_mask], indices[span_mask]


def distance2segments_2(neut_pos, segments, cst, span):
    """create inverse segments and recursively search neutral sites in them"""
    # load genetic map class
    gmp = GeneticMap(cst.chrom, cst.gmap_files)

    # create inner segments
    starts = np.concatenate(([0], segments[:, 1]))
    ends = np.concatenate((segments[:, 0], [cst.chlen]))
    segs = np.column_stack((starts, ends))

    # find which inner segment neutral site lands in
    indices, distance = [], []
    for (i, npos) in enumerate(neut_pos):
        istart, iend = recursive_search(npos, segs)
        gstart = gmp.interp_gpos(istart)
        gend = gmp.interp_gpos(iend)
        gpos = gmp.interp_gpos(npos)
        dist5 = gpos - gend
        dist3 = gpos - gstart
        assert dist5 <= 0 <= dist3
        if abs(dist5) > dist3:
            distance.append(dist3)
        else:
            distance.append(dist5)
        indices.append(i)

    indices = np.array(indices)
    distance = np.array(distance)
    span_mask = (abs(distance) <= span)

    return distance[span_mask], indices[span_mask]


def bin_distances(cst, chrom_data):
    """
    Bin data by distance to nearest focal site
    :param cst: DataStruct
    :param chrom_data: [distance, div, pi, pred] array for a single chromosome
    :return binned_data: [div, pi, pred, count] summed into bins defined by {width, span} params from cst,
    where count is the number of sites in each of the bin sums 
    """
    # get the span of the collated plot and the width of bin windows from cst
    span, width = cst.vars.collated_span, cst.vars.collated_bin

    # bin the distances for each neutral site to collect pi and pred values
    cnt, bins = np.histogram(a=chrom_data[0], bins=np.arange(-span, span+width, width))

    # store pi, div, pred, count for each bin in ar
    binned_data = np.zeros(shape=[len(cnt), 4])

    i = 0
    for idx in xrange(len(cnt)):

        # define upper data index by the count in the bin
        j = i + cnt[idx]

        # if theres data, take sum of [div, pi, pred, count] for current window
        if cnt[idx]:
            binned_data[idx, :3] = np.sum(chrom_data[1:, i:j], axis=1)
            binned_data[idx, 3] = cnt[idx]
        # empty bins get 0s
        else:
            binned_data[idx] = [0.0, 0.0, 0.0, 0.0]

        # set lower range to current upper range
        i = j

    return binned_data


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
