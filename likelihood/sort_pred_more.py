__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from precalc.lh_inputs import summary_stats
from classes.knowngene import complement_strand
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import load_saved, adjust_arrays
from classes.runstruct import ChromStruct, root_dir, FileStruct, \
    human_autosomes, cst_from_fldr


def basic_sort(nt, ns, nu, gc, cn, cm, il, pred, num):
    """Sort data arrays by predictions and average results into bins"""
    # get sorting indices based on prediction value
    sidx = np.argsort(pred)
    # sort predictions and div/poly data
    all_arrs = nt, ns, nu, gc, cn, cm, il, pred
    nt, ns, nu, gc, cn, cm, il, pred = [a[sidx].astype('f8') for a in all_arrs]

    # calculate mean div
    meandiv = np.average(nu, weights=ns)

    # get start and end indices per partition
    idx = sortbin_edges(sites=ns, numbins=num)

    # gather data summaries
    sorted_array = []
    for (i, j) in idx:
        # calculate each statistic for the present bin
        norm_div = np.average(nu[i:j], weights=ns[i:j]) / meandiv
        scaled_pi = np.sum(nt[i:j, 1]) / (np.sum(nt[i:j]) * norm_div)
        # use weighted mean for pred
        prediction = np.average(pred[i:j], weights=ns[i:j])
        # GC is already weighted per neutral site
        mean_gc = np.mean(gc[i:j])
        # conservation and cM/Mb are not weighted for sites
        mean_cn = np.average(cn[i:j], weights=ns[i:j])
        mean_cm = np.average(cm[i:j], weights=ns[i:j])
        mean_il = np.average(il[i:j], weights=ns[i:j])

        sorted_array.append([norm_div, scaled_pi, prediction,
                             mean_gc, mean_cn, mean_cm, mean_il])

    return np.array(sorted_array)


def archaic_sort(ai, ns, pred, num):
    """Sort data arrays by predictions and average results into bins"""
    # get sorting indices based on prediction value
    sidx = np.argsort(pred)
    # sort predictions and div/poly data
    all_arrs = ai, ns, pred
    ai, ns, pred = [a[sidx].astype('f8') for a in all_arrs]

    # get start and end indices per partition
    idx = sortbin_edges(sites=ns, numbins=num)

    # calculate mean introgression levels in each bin weighted by sites
    sorted_ai = []
    for (i, j) in idx:
        prediction = np.average(pred[i:j], weights=ns[i:j])
        # arch = np.average(ai[i:j], weights=ns[i:j])
        arch = np.sum(ai[i:j]) / np.sum(ns[i:j])
        sorted_ai.append((prediction, arch))

    return np.array(sorted_ai)


def segment_sort(cst, ns, sg, pred, num=100):
    """Sort segment arrays by prediction and save one for each bin"""
    # create new save dir if needed
    bin_dir = root_dir + '/compress/sorted_bins'
    if not os.path.isdir(bin_dir):
        os.mkdir(bin_dir)

    # get sorting indices based on prediction value
    sidx = np.argsort(pred)

    # sort segments and number of sites per segment
    sg, ns = sg[sidx], ns[sidx]

    # get start and end indices per partition
    idx = sortbin_edges(sites=ns, numbins=num)

    # for each percentile, create an array of the included sites
    bin_num = 0
    for (i, j) in idx:
        seg_group = sg[i:j]
        seg_sites = ns[i:j]
        # f_save = bin_dir + '/sorted_bin{}.npz'.format(bin_num)
        f_save = bin_dir + '/sorted_bin{}.{}.npz'.format(bin_num, cst.tkn)
        np.savez_compressed(f_save, sg=np.column_stack((seg_group, seg_sites)))
        # np.savez_compressed(f_save, sg=seg_group)
        bin_num += 1

    return None


def snptype_sort(snp_dict, ns, pred, num, notsorted=False):
    """Sort segment arrays by prediction and save one for each bin"""
    # get sorting indices based on prediction value
    sidx = np.argsort(pred)
    # sort segments and number of sites per segment
    ns = ns[sidx]

    # for each SNP type, sort the poly data
    snp_types = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT']
    labels = []
    for snptype in snp_types:
        # get complement SNP type
        complement = complement_strand(snptype)
        # create label for current SNP type
        lbl = '{}_{}'.format(snptype, complement)
        labels.append(lbl)
        if notsorted:
            snp_dict[lbl] = snp_dict[lbl][sidx]

    # get start and end indices per partition
    idx = sortbin_edges(sites=ns, numbins=num)

    # for each percentile, get mean pi per SNP type
    sort_array = []
    for (i, j) in idx:
        row = []
        for lbl in labels:
            nt = snp_dict[lbl]
            pi = np.sum(nt[i:j, 1]) / np.sum(nt[i:j])
            row.append(pi)
            # add number of sites for index
        num_sites = ns[i:j].sum()
        num_segs = j-i
        row.append(num_sites)
        row.append(num_segs)
        sort_array.append(row)

    sort_array = np.array(sort_array)

    return sort_array


def load_more(cst):
    assert isinstance(cst, ChromStruct)
    gc, cn, cm, il = [], [], [], []
    for ch in cst.chroms:
        cst.chrom = ch
        gc.append(np.load(cst.gc_files)['gc'])
        cn.append(np.load(cst.cn_files)['cn'])
        cm.append(np.load(cst.cm_files)['cm'])
        il.append(np.load(cst.il_files)['il'])
    gc = np.concatenate(gc)
    cn = np.concatenate(cn)
    cm = np.concatenate(cm)
    il = np.concatenate(il)

    return gc, cn, cm, il


def load_seg(cst):
    """load segment data for each chrom"""
    sg_dir = root_dir + '/compress/sg_arr/'
    sg = []
    for ch in human_autosomes:
        f_load = sg_dir + '{}.sg_arr.{}.npz'.format(ch, cst.tkn)
        sg.append(np.load(f_load)['sg'])
    sg = np.concatenate(sg)

    return sg


def load_snptypes(cst, msk):
    """load snp type data for each snp type"""
    # create save directory for SNP types
    snp_dir = root_dir + '/compress/snp_type/'
    # create format for saved files
    # load_fmt = snp_dir + '{}.{}.anc_corrected_3.npz'
    load_fmt = snp_dir + '{}.{}.{}.fix_anc.npz'

    # for each SNP type, get compressed poly data
    snp_dict = {}
    snp_types = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT']
    for snptype in snp_types:
        # get complement SNP type
        complement = complement_strand(snptype)
        # create label for current SNP type
        lbl = '{}_{}'.format(snptype, complement)
        # load all data for each snp type
        cur_data = []
        for ch in human_autosomes:
            f_load = load_fmt.format(ch, cst.tkn, lbl)
            data = np.load(f_load)['arr_0']
            cur_data.append(data)
        snp_dict[lbl] = np.concatenate(cur_data)[msk]

    return snp_dict


def load_archaic(cst, pop_id):
    ai = []
    for ch in human_autosomes:
        cst.chrom = ch
        f_load = cst.ai_files.replace('ai.', 'ai.{}.'.format(pop_id))
        ai.append(np.load(f_load)['ai'])
    ai = np.concatenate(ai)

    return ai


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


def main_1():
    if len(argv) != 2:
        print 'usage: sort_pred_more <folder_name>'
        exit(1)

    fldr = argv[1]
    fdir = root_dir + '/result/final_files/{}'.format(fldr)
    cst = cst_from_fldr(fldr)
    cst.files = FileStruct(cst.dict)
    # tau = cst.params[-1]

    # fldr = 'cadd93_align8'
    # fdir = root_dir + '/result/final_files/{}'.format(fldr)
    # # cst = cst_from_fldr(fldr)
    # # cst.files = FileStruct(cst.dict)
    # # fldr = argv[1]
    # # cst = cst_from_fldr(fldr)
    # # cst.files = FileStruct(cst.dict)
    # f_init = argv[1]
    # cst = ChromStruct('chr1', init=f_init)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # load more data
    gc, cn, cm, il = load_more(cst)

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
    msk &= mnu1[:,0]

    gc, cn, cm, il = gc[msk], cn[msk], cm[msk], il[msk]

    # use constant mutation rate for predictions!
    if fldr == 'cadd93_extel_rm_CG':
        nu_pred = np.full(len(nu), np.average(nu, weights=dv, axis=0))
        use_params = cst_from_fldr('cadd93_bth600').params
        print nu_pred.mean()
        print use_params
    else:
        nu_pred = np.full(len(nu), cst.stat.meandiv)
        use_params = cst.params

    # # use standard params for BGC or CpG
    # if ('BGC' in fldr) or ('CpG' in fldr):
    #     print 'Using default params for predictions'
    #     use_params = cst_from_fldr('cadd93_new_nu').params
    #     sort_file = fdir + '/predsort_more_stdsort.txt'
    # else:
    #     use_params = cst.params
    #     sort_file = fdir + '/predsort_more.txt'
    # use_params[-1] = tau

    # set manually for new version
    # get params from cadd93 run
    # use_params = cst_from_fldr('cadd93').params
    # use_params = cst.params
    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(use_params, cst, nu_pred, bs, cs)
    print pred.mean()
    # get sorted array, default n=100
    # n = 100
    for n in [100, 2000]:
    # n = 100
        sort_file = fdir + '/sort_gc_cm_cn_il_n{}.txt'.format(n)
        sarr = basic_sort(nt, dv, nu, gc, cn, cm, il, pred, n)
        np.savetxt(sort_file, sarr)


def main_2():
    """CHROMOSOME SEGMENT ARRAY SORTING"""
    if len(argv) != 2:
        print 'usage: sort_pred_more <folder_name>'
        exit(1)

    fldr = argv[1]
    cst = cst_from_fldr(fldr)
    cst.files = FileStruct(cst.dict)
    # f_init = argv[1]
    # cst = ChromStruct('chr1', init=f_init)

    # load complete array data
    _, bs, cs, nu, nt, dv, pl = load_saved(cst)

    # load more data
    sg = load_seg(cst)
    # print sg.shape

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    cutoff = 0.25
    mnu1 = (nu <= cutoff)

    # mask and rescale maps
    bs, cs, nu, nt, dv, pl, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)
    # print msk.shape
    # print msk.sum()

    # TEMPORARY FIX FOR ABERRANT SUBSTITUTION RATE ESTIMATES!!!
    mnu2 = (nu <= cutoff)
    if bs is not None:
        bs = bs[mnu2]
    if cs is not None:
        cs = cs[mnu2]
    nu, nt, dv, pl = [a[mnu2] for a in nu, nt, dv, pl]
    msk &= mnu1[:,0]

    sg = sg[msk]

    # use constant mutation rate for predictions!
    if fldr == 'cadd93_extel_rm_CG':
        nu_pred = np.full(len(nu), np.average(nu, weights=dv, axis=0))
        use_params = cst_from_fldr('cadd93_bth600').params
    else:
        nu_pred = np.full(len(nu), cst.stat.meandiv)
        use_params = cst.params

    # get params from cadd93 run
    # use_params = cst_from_fldr('cadd93').params
    # pred = predicted_pi(use_params, cst, nu_pred, bs, cs)

    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(use_params, cst, nu_pred, bs, cs)

    # print 'mean pred = {}'.format(pred.mean())

    # get sorted array, default n=100
    n = 2000
    # n = 100
    segment_sort(cst, dv, sg, pred, num=n)


def main_3():
    """SNP TYPE SORTING"""
    if len(argv) != 2:
        print 'usage: sort_pred_more <folder_name>'
        exit(1)

    # fldr = 'cadd93_align8'
    # fdir = root_dir + '/result/final_files/{}'.format(fldr)
    # cst = cst_from_fldr(fldr)
    # cst.files = FileStruct(cst.dict)

    # f_init = argv[1]
    # cst = ChromStruct('chr1', init=f_init)

    fldr = argv[1]
    cst = cst_from_fldr(fldr)
    cst.files = FileStruct(cst.dict)
    fdir = root_dir + '/result/final_files/{}'.format(fldr)

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

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
    msk &= mnu1[:,0]

    # load SNP dict
    snp_dict = load_snptypes(cst, msk)

    # use constant mutation rate for predictions!
    if fldr == 'cadd93_extel_rm_CG':
        nu_pred = np.full(len(nu), np.average(nu, weights=dv, axis=0))
        use_params = cst_from_fldr('cadd93_bth600').params
    else:
        nu_pred = np.full(len(nu), cst.stat.meandiv)
        use_params = cst.params
    # get params from cadd93 run
    # use_params = cst_from_fldr('cadd93').params
    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(use_params, cst, nu_pred, bs, cs)

    # get sorted array, default n=100
    # n = 100
    sort_bins = [100, 2000]
    for (i, n) in enumerate(sort_bins):
        if i == 0:
            sarr = snptype_sort(snp_dict, dv, pred, n, notsorted=True)
        else:
            sarr = snptype_sort(snp_dict, dv, pred, n, notsorted=False)
        sort_file = fdir + '/sort_snptype_fix_anc_n{}.txt'.format(n)
        np.savetxt(sort_file, sarr)


def main_4():
    """ARCHAIC INTROGRESSION SORTING"""
    if len(argv) != 2:
        print 'usage: sort_pred_more <folder_name>'
        exit(1)
    
    # load cst from save folder
    fldr = argv[1]
    fdir = root_dir + '/result/final_files/{}'.format(fldr)
    cst = cst_from_fldr(fldr)
    cst.files = FileStruct(cst.dict)  # update for added FileStruct variables

    # load complete array data
    sg, bs, cs, nu, nt, dv, pl = load_saved(cst)

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
    msk &= mnu1[:,0]

    # use constant mutation rate for predictions!
    if fldr == 'cadd93_extel_rm_CG':
        nu_pred = np.full(len(nu), np.average(nu, weights=dv, axis=0))
        use_params = cst_from_fldr('cadd93_bth600').params
    else:
        nu_pred = np.full(len(nu), cst.stat.meandiv)
        use_params = cst.params

    # calculate predicted pi from arrays and average of top 3 params
    pred = predicted_pi(use_params, cst, nu_pred, bs, cs)

    # get sorted array, default n=100
    # n = 100
    # pop_id = 'CHBS'
    for pop_id in ['CEU', 'CHBS']:
        # load archaic introgression and mask
        ai = load_archaic(cst, pop_id)[msk]
        for n in [100, 2000]:
            # n = 100
            sarr = archaic_sort(ai, dv, pred, n)
            f_save = fdir + '/predsort_archaic_{}_n{}.txt'.format(pop_id, n)
            np.savetxt(f_save, sarr)


if __name__ == '__main__':
    # main_1()  # sort extra data: cM, cons, GC-content, CpG island content
    main_2()  # create B bins
    # main_3()  # sort SNP types
    # main_4()  # sort archaic introgression
