__author__ = 'davidmurphy'

import os
import numpy as np
import scipy.stats as sts
from sys import argv
from itertools import izip
from precalc.lh_inputs import load_saved
from data_processing.functions import calc_pi
from data_processing.data_tools import snpcount
from data_processing.precalc_tools import neutral_snps
from classes.runstruct import RunStruct, ChromStruct, root_dir
from scipy.stats import linregress
from scipy.stats import t as ttest


def re_save(npy_file):
    a = np.load(npy_file)
    fz = npy_file.replace('.npy', '.npz')
    if 'pivals' in npy_file:
        np.savez_compressed(fz, pi=a)
    if 'pdvals' in npy_file:
        np.savez_compressed(fz, pd=a)
    if 'nvals' in npy_file:
        np.savez_compressed(fz, n=a)
    if 'bvals' in npy_file:
        np.savez_compressed(fz, b=a)


def get_pi_vals(cst, nmsk, pop):
    """get neutral positions and calculate pi for those with SNPs"""
    # get neutral SNP file for the correct pop
    fsnp = cst.snp_files.replace(cst.neut, pop)
    print fsnp

    # get neutral SNP data
    pos, p, q = neutral_snps(nmsk, np.column_stack(snpcount(fsnp))).T

    # calculate pi
    pi = calc_pi(p + q, q)

    return pos, pi


def get_pop_div(cst, nmsk, pop_1, pop_2):
    """get population divergence at neutral sites with one or more SNPs"""
    assert isinstance(cst, ChromStruct)

    # set SNP file paths
    f_1 = cst.snp_files.replace(cst.neut, pop_1)
    f_2 = cst.snp_files.replace(cst.neut, pop_2)

    # get neutral snps
    s_1 = neutral_snps(nmsk, np.column_stack(snpcount(f_1)))
    s_2 = neutral_snps(nmsk, np.column_stack(snpcount(f_2)))

    # fill spaces in each SNP set with monomorphic sites
    # s_1, s_2 = join_pop_snps(s_1, s_2)
    assert np.all(s_1[:, 0] == s_2[:, 0])

    # calculate the het pairs and hom pairs between the populations
    dhet = (s_1[:, 1] * s_2[:, 2]) + (s_1[:, 2] * s_2[:, 1])
    dhom = (s_1[:, 1] * s_2[:, 1]) + (s_1[:, 2] * s_2[:, 2])

    # calculate the pop divergence
    dpop = 1.0 * dhet / (dhom + dhet)

    # return SNP positions and pop divergence
    return s_1[:, 0], dpop


def join_pop_snps(snp_1, snp_2):
    # combine all positions into single array
    all_pos = np.unique(np.concatenate((snp_1[:, 0], snp_2[:, 0])))
    print all_pos.size, snp_1[:, 0], snp_2[:, 0]

    # get the sample size for filling in gaps
    sample_1 = np.sum(snp_1[0, 1:])
    sample_2 = np.sum(snp_2[0, 1:])
    # print sample_1, sample_2

    # get indices of existing SNPs for each pop
    si_1 = np.in1d(all_pos, snp_1[:, 0])
    si_2 = np.in1d(all_pos, snp_2[:, 0])

    # refill pop 1 arrays with existing poly and mono sites
    s_1 = np.zeros(shape=(all_pos.size, 3))
    # fill in all positions
    s_1[:, 0] = all_pos
    # fill in existing SNPs
    s_1[:, 1:][si_1] = snp_1[:, 1:]
    # fill in monomorphic sites with sample size
    s_1[:, 1:][~si_1] = [sample_1, 0]

    # refill pop 2 arrays with existing poly and mono sites
    s_2 = np.zeros(shape=(all_pos.size, 3))
    # fill in all positions
    s_2[:, 0] = all_pos
    # fill in existing SNPs
    s_2[:, 1:][si_2] = snp_2[:, 1:]
    # fill in monomorphic sites with sample size
    s_2[:, 1:][~si_2] = [sample_2, 0]

    return s_1, s_2


def bin_sorted_vals(pop1, pop2, nbins=100):
    # set file path
    rdir = 'ape95_minbs_01'
    fdir = root_dir + '/result/final_files/{}/popdiv'.format(rdir)

    # set file names
    fb_arra = '{}/bvals.npz'.format(fdir)
    p1_arra = '{}/{}.pivals.npz'.format(fdir, pop1)
    p2_arra = '{}/{}.pivals.npz'.format(fdir, pop2)
    pd_arra = '{}/{}.{}.pdvals.npz'.format(fdir, pop1, pop2)
    nn_arra = '{}/nvals.npz'.format(fdir)

    # re-save the remaining .npy versions of files
    for f in [p1_arra, p2_arra, pd_arra]:
        if not os.path.isfile(f):
            re_save(f.replace('.npz', '.npy'))

    # load files
    bv = np.load(fb_arra)['b']
    pi1 = np.load(p1_arra)['pi']
    pi2 = np.load(p2_arra)['pi']
    pd = np.load(pd_arra)['pd']
    nn = np.load(nn_arra)['n']

    # set step size by number of bins
    step = len(bv) / nbins

    # sort data by b values
    si = np.argsort(bv)
    bv, pi1, pi2, pd, nn = bv[si], pi1[si], pi2[si], pd[si], nn[si]

    # save mean data values in bins
    bmean, pimean1, pimean2, pdmean = [], [], [], []
    for i in range(0, bv.size, step):
        j = i + step
        bmean.append(np.average(bv[i:j], weights=nn[i:j]))
        pimean1.append(np.average(pi1[i:j], weights=nn[i:j]))
        pimean2.append(np.average(pi2[i:j], weights=nn[i:j]))
        pdmean.append(np.average(pd[i:j], weights=nn[i:j]))

    # convert results to arrays
    pimean1 = np.array(pimean1)
    pimean2 = np.array(pimean2)
    pdmean = np.array(pdmean)
    bmean = np.array(bmean)

    # create paths to save files
    bsave = '{}/n{}.bmean.txt'.format(fdir, nbins)
    pisave1 = '{}/{}.n{}.pimean.txt'.format(fdir, pop1, nbins)
    pisave2 = '{}/{}.n{}.pimean.txt'.format(fdir, pop2, nbins)
    pdsave = '{}/{}.{}.n{}.pdmean.txt'.format(fdir, pop1, pop2, nbins)

    # save each file if it does not exist
    if not os.path.isfile(bsave):
        np.savetxt(bsave, bmean)
    if not os.path.isfile(pisave1):
        np.savetxt(pisave1, pimean1)
    if not os.path.isfile(pisave2):
        np.savetxt(pisave2, pimean2)
    if not os.path.isfile(pdsave):
        np.savetxt(pdsave, pdmean)


def sort_pop_data(pop_1, pop_2):
    # use standard run as default
    fdir = root_dir + '/result/final_files/ape95_minbs_01/'
    flist = [f for f in os.listdir(fdir) if f.endswith('composite.txt')]
    assert len(flist) == 1

    # create paths to rsq and final init file
    f_init = fdir + flist[0]
    cst = ChromStruct(chrom='chr1', init=f_init)

    # get data from each chrom
    b_vals, pd_vals, p1_vals, p2_vals, n_vals = [], [], [], [], []
    for ch in cst.chroms:
        # reset chrom
        cst.chrom = ch

        # load array data for chrom
        sg, bs, _, nu, _, _, _ = load_saved(cst, chroms=[ch])
        bs *= np.log1p(-1.0 / cst.bscl)
        bs = np.maximum(np.log(0.01), bs)
        pos = np.cumsum(sg)

        # load neutral mask
        nmsk = np.load(cst.neut_masks)['neutmask']

        # get pop div for two pops
        spos, pdiv = get_pop_div(cst, nmsk, pop_1, pop_2)

        # get pi values for each pop
        sp_1, pi_1 = get_pi_vals(cst, nmsk, pop_1)
        sp_2, pi_2 = get_pi_vals(cst, nmsk, pop_2)

        # check that all positions are the same
        print spos.size, sp_1.size, sp_2.size
        assert np.all(spos == sp_1) and np.all(spos == sp_2)

        # get positions of neutral sites (0-based)
        nidx = np.where(nmsk != 0)[0]

        # create new arrays to store pi and pop div values
        popdiv = np.zeros(shape=nidx.size)
        mpi_1 = np.zeros(shape=nidx.size)
        mpi_2 = np.zeros(shape=nidx.size)

        # find position with polymorphism
        dmsk = np.in1d(nidx, spos)
        assert np.sum(dmsk) == len(spos)

        # fill polymorphism or divergence data at poly sites
        popdiv[dmsk] = pdiv
        mpi_1[dmsk] = pi_1
        mpi_2[dmsk] = pi_2

        # get indices of segments where neutral polymorphic hits occurred
        sidx = np.concatenate(([0], np.searchsorted(nidx, pos)))

        # create new lists to store the number of sites and mean values per seg
        num, pop_div_mean, pi_1_mean, pi_2_mean = [], [], [], []
        for (i, j) in izip(sidx[:-1], sidx[1:]):
            n = j - i
            if n > 0:
                pd_mean = np.mean(popdiv[i:j])
                p1_mean = np.mean(mpi_1[i:j])
                p2_mean = np.mean(mpi_2[i:j])
            else:
                pd_mean = 0
                p1_mean = 0
                p2_mean = 0
            num.append(n)
            pop_div_mean.append(pd_mean)
            pi_1_mean.append(p1_mean)
            pi_2_mean.append(p2_mean)

        # convert lists of mean values to arrays
        num = np.array(num)
        pop_div_mean = np.array(pop_div_mean)
        pi_1_mean = np.array(pi_1_mean)
        pi_2_mean = np.array(pi_2_mean)

        # get the bmap for each chrom
        params = cst.params
        uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
        bwt = uvec / cst.fixed.u_fix
        bsx = np.exp(np.dot(bs, bwt))
        assert bsx.size == num.size == pop_div_mean.size == pi_1_mean.size == \
               pi_2_mean.size

        # adjust pi values by divergence estimate
        pop_div_mean /= (nu[:, 0] / 0.1356351719743293)
        pi_1_mean /= (nu[:, 0] / 0.1356351719743293)
        pi_2_mean /= (nu[:, 0] / 0.1356351719743293)

        # append values where number of sites > 0
        b_vals.append(bsx[num > 0])
        pd_vals.append(pop_div_mean[num > 0])
        p1_vals.append(pi_1_mean[num > 0])
        p2_vals.append(pi_2_mean[num > 0])
        n_vals.append(num[num > 0])

    # concatenate all data across chroms
    b_vals = np.concatenate(b_vals)
    pd_vals = np.concatenate(pd_vals)
    p1_vals = np.concatenate(p1_vals)
    p2_vals = np.concatenate(p2_vals)
    n_vals = np.concatenate(n_vals)

    # save pop div file by default
    sdir = fdir + 'popdiv/'
    pop_div_file = '{}{}.{}.pdvals.npz'.format(sdir, pop_1, pop_2)
    np.savez_compressed(pop_div_file, pd=pd_vals)
    # only save the pi files if they are not duplicates
    pi1_save = sdir + pop_1 + '.pivals.npz'
    pi2_save = sdir + pop_2 + '.pivals.npz'
    if not os.path.isfile(pi1_save):
        np.savez_compressed(pi1_save, pi=p1_vals)
    if not os.path.isfile(pi2_save):
        np.savez_compressed(pi2_save, pi=p2_vals)
    # only need to save num file once
    nn_save = sdir + 'nvals.npz'
    if not os.path.isfile(nn_save):
        np.savez_compressed(nn_save, n=n_vals)
    # only need to save b files once
    bb_save = sdir + 'bvals.npz'
    if not os.path.isfile(bb_save):
        np.savez_compressed(bb_save, b=b_vals)


def compare_pop_data(pop_1, pop_2, nbins=1000):
    # set file path
    rdir = 'ape95_minbs_01'
    fdir = root_dir + '/result/final_files/{}/popdiv'.format(rdir)

    # create paths to binned data files
    pisave1 = '{}/{}.n{}.pimean.txt'.format(fdir, pop_1, nbins)
    pisave2 = '{}/{}.n{}.pimean.txt'.format(fdir, pop_2, nbins)
    pdsave = '{}/{}.{}.n{}.pdmean.txt'.format(fdir, pop_1, pop_2, nbins)
    f_bval = '{}/n{}.bmean.txt'.format(fdir, nbins)

    # subtract div-pi intercepts for middle 50% of bins
    i = int(nbins * 0.05)
    j = int(nbins * 0.95)
    # n = j - i

    # load data
    p1 = np.loadtxt(pisave1)[i:j]
    b = np.loadtxt(f_bval)[i:j]

    # output format
    # fmt = '{:.0f},{:.0f},{:.0f}'
    fmt = '{:.0f} [{:.0f},{:.0f}]'

    if pop_1 == pop_2:
        res1 = get_ci_yr(b, p1)
        # res1 = [r * 30.0 / (2 * 1.2e-8) for r in mean1, lo1, hi1]
        output1 = fmt.format(*res1)

        return output1

    else:
        fmt = '{:.0f}'
        # loaded additional data
        p2 = np.loadtxt(pisave2)[i:j]
        pd = np.loadtxt(pdsave)[i:j]
        b = np.loadtxt(f_bval)[i:j]

        # subtract population 1 and then population 2
        sub1 = pd - p1
        sub2 = pd - p2
        # use intercept formula on subtraction
        # res1 = get_ci_yr(b, sub1)
        res1 = get_ci_yr(b, pd)[0] - get_ci_yr(b, p1)[0]

        # get mean, stderr and 95% ci from each subtraction
        # mean1 = np.mean(sub1)
        # serr1 = np.std(sub1)
        # ci1 = serr1 * sts.t.ppf(0.975, n-1)
        # lo1, hi1 = mean1-ci1, mean1+ci1
        # put in units of years split

        # use intercept formula on subtraction
        # mean2 = np.mean(sub2)
        # serr2 = np.std(sub2)
        # ci2 = serr2 * sts.t.ppf(0.975, n-1)
        # lo2, hi2 = mean2-ci2, mean2+ci2
        # res2 = get_ci_yr(b, sub2)
        res2 = get_ci_yr(b, pd)[0] - get_ci_yr(b, p2)[0]

        # res2 = [r * 30.0 / (2 * 1.2e-8) for r in mean2, lo2, hi2]

        output1 = fmt.format(res1)
        output2 = fmt.format(res2)

        return output1, output2


def ci_intercept(b_0, k, n, se, alpha):
    """
    https://www.tutorialspoint.com/statistics
    /regression_intercept_confidence_interval.htm
    """
    t1 = 1 - ((1-alpha) / 2.0)
    t2 = n - k - 1
    r_high = b_0 + ttest.ppf(t1, t2) * se
    r_low = b_0 - ttest.ppf(t1, t2) * se

    return r_low, r_high


def get_intercept_ci(xi, yi):
    n = len(xi)
    slope, intercept, rvalue, pvalue, stderr = linregress(xi, yi)
    lfunc = lambda x: slope * x + intercept
    yi_err = yi - lfunc(xi)
    sigma_sq = np.sum(yi_err**2) / (n-2)
    xi_mean = np.mean(xi)
    xi_var = np.sum((xi-xi_mean)**2)
    b0_se = np.sqrt(sigma_sq * (1.0/n + xi_mean**2 / xi_var))
    ilow, ihigh = ci_intercept(intercept, 1, n, b0_se, 0.95)
    return intercept, ilow, ihigh


def get_ci_yr(xi, yi):
    ci_data = get_intercept_ci(xi, yi)
    yr_int, yr_low, yr_high = [x * 30.0 / (2*1.2e-8) for x in ci_data]
    # return yr_low, yr_int, yr_high
    return yr_int, yr_low, yr_high


def make_table_1():
    rdir = 'ape95_minbs_01'
    fdir = root_dir + '/result/final_files/{}/popdiv'.format(rdir)
    f_pops = '{}/super_pops.txt'.format(fdir)
    f_tabl = '{}/table_1.txt'.format(fdir)
    # f_tabl = '{}/table_1_viewable.txt'.format(fdir)

    # get population labels
    pops = []
    sups = []
    with open(f_pops, 'r') as f:
        for pop in f:
            pops.append(pop.split()[0])
            sups.append(pop.split()[1])

    pop_sup = ['{}-{}'.format(p, s) for p, s in zip(pops, sups)]
    with open(f_tabl, 'w') as f:
        # create header row with pop labels
        f.write('X\t{}\n'.format('\t'.join(pop_sup)))

        # go through each population comparison and save to table
        for i in range(26):
            row = []
            pop1 = pops[i]
            for j in range(26):
                pop2 = pops[j]
                if pop1 == pop2:
                    output = compare_pop_data(pop1, pop2)
                else:
                    try:
                        output = compare_pop_data(pop1, pop2)[0]
                    except IOError:
                        output = compare_pop_data(pop2, pop1)[1]
                row.append(output)

            # write row of results
            f.write('{}\t{}\n'.format(pop_sup[i], '\t'.join(row)))


def make_table_2():
    rdir = 'ape95_minbs_01'
    fdir = root_dir + '/result/final_files/{}/popdiv'.format(rdir)
    f_pops = '{}/super_pops.txt'.format(fdir)
    f_tabl = '{}/table_2.txt'.format(fdir)

    # get population labels
    pops = []
    sups = []
    with open(f_pops, 'r') as f:
        for pop in f:
            pops.append(pop.split()[0])
            sups.append(pop.split()[1])

    # get start/end indices for middle 50% of values
    nb = 1000
    ii = int(nb * 0.05)
    jj = int(nb * 0.95)

    # load b values
    f_bval = '{}/n{}.bmean.txt'.format(fdir, nb)
    b = np.loadtxt(f_bval)[ii:jj]

    # start table file
    pop_sup = ['{}-{}'.format(p, s) for p, s in zip(pops, sups)]
    with open(f_tabl, 'w') as f:
        # create header row with pop labels
        f.write('X\t{}\n'.format('\t'.join(pop_sup)))

        # go through each population comparison and save to table
        for i in range(len(pops)-1):
            pop1 = pops[i]
            # pad row with N/A for redundant comparisons
            row = ['X'] * (i+1)
            for j in range(i+1, len(pops), 1):
                pop2 = pops[j]
                fmt = '{}/{}.{}.n{}.pdmean.txt'
                try:
                    f_dval = fmt.format(fdir, pop1, pop2, nb)
                    pd = np.loadtxt(f_dval)[ii:jj]
                except IOError:
                    f_dval = fmt.format(fdir, pop2, pop1, nb)
                    pd = np.loadtxt(f_dval)[ii:jj]

                intc, lo, hi = get_ci_yr(b, pd)
                row.append('{:.0f} [{:.0f},{:.0f}]'.format(intc, lo, hi))

            # write row of results
            f.write('{}\t{}\n'.format(pop_sup[i], '\t'.join(row)))


def main():
    make_table_1()
    make_table_2()
    # if len(argv) != 3:
    #     print 'usage: pop_div_sort <pop_1> <pop_2>'
    #     exit(1)
    #
    # # get populations form command line
    # pop_1, pop_2 = argv[1:]
    # # sort_pop_data(pop_1, pop_2)
    # # bin_sorted_vals(pop_1, pop_2, nbins=1000)
    # compare_pop_data(pop_1, pop_2)


if __name__ == '__main__':
    main()
