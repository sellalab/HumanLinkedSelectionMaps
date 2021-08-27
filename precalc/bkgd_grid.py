import numpy as np
from classes.runstruct import ChromStruct
from classes.geneticmap import GeneticMap
from classes.annosegments import AnnoSegments
from classes.bkgdcalculator import BkgdCalculator
from data_processing.data_tools import randint_unique

# import matplotlib.pyplot as plt
# import seaborn


__author__ = 'davidmurphy'


def bkgd_grid(cst, anno, coef, err):
    """get maximally spaced points where b changes by less than err"""
    # use BkgdCalculator to summing derivatives at each step
    assert cst.chrom == anno.chrom
    udel = cst.fixed.u_fix
    bkc = BkgdCalculator(coef, udel, anno)

    # get the lower and upper bound of the chromosome
    y_min = bkc.gmap.interp_gpos(1.0)
    y_max = bkc.gmap.interp_gpos(cst.chlen)

    # initial step d_0 = 1/5 the gdist to first cons segment
    # d_0 = 0.2 * (anno.rstart.min() - y_min)
    d_0 = 1e-3

    # start the algorithm at the beginning of the chromosome
    y = y_min
    grid_points = []
    while y < y_max:
        step = minimum_stepsize(bkc, y, d_0, err)
        y += step
        grid_points.append(y)

    return np.array(grid_points)


def minimum_stepsize(bkc, y, d_0, err):
    # select points behind y, sum 1st derivative terms
    ma = (bkc.gcons < y)
    x_a = y - bkc.gcons[ma]
    db_a = bkc.bkgd_sum_der1(0.0, x_a)

    # select points beyond y+d_0, sum 1st derivative terms
    mb = (bkc.gcons > (y + d_0))
    x_b = bkc.gcons[mb] - y - d_0
    db_b = bkc.bkgd_sum_der1(0.0, x_b)

    # select points in the interval [y, y+d_0], sum 1st derivative terms
    m_ab = ~(ma | mb)
    # use WORST case scenario for middle sites and assume y=x
    # x_ab = np.zeros(shape=np.sum(m_ab))
    # use MIDPOINT approximation for inner sites
    x_ab = abs(bkc.gcons[m_ab] - y - 0.5*d_0)
    db_ab = bkc.bkgd_sum_der1(0.0, x_ab)

    # calculate upper bound of the 1st derivative in [y, y+d_0]
    max_db = db_a + db_b + db_ab
    # maximum step length: err / (db/dy)
    min_dy = err / max_db

    # calculate genetic map position after moving 1bp for lower bound
    dy_bp = bkc.gmap.interp_gpos(bkc.gmap.interp_pos(y)+1) - y

    # bound the step between 1bp and d_0
    step = min(d_0, max(dy_bp, min_dy))

    return step


def bterm_contribution(cst, anno, coef):
    """break down the contribution of each conserved site to the bsum/dbsum"""
    # use BkgdCalculator to summing derivatives at each step
    assert cst.chrom == anno.chrom
    udel = cst.fixed.u_fix
    bkc = BkgdCalculator(coef, udel, anno)

    # take some random conserved positions
    n = anno.num_bases
    r_idx = randint_unique(5e3, n)
    # r_idx = randint_unique(1e3, start=int(0.25 * n), end=int(0.75 * n))
    r_gpos = bkc.gcons[r_idx]

    # get b terms as a vector
    fraction = []
    distance = []
    bp_cm = []
    density = []
    for gy in r_gpos:
        bv = bkc.bkgd_sum_chrom(gy, bsum=False)
        # bv = bkc.bkgd_sum_chrom_der1(gy, dbsum=False)
        # cumulative sum of b/db terms
        fcum = np.cumsum(bv / np.sum(bv))
        # get mask for sites having middle 98% of effect
        mi = (fcum > 0.01) & (fcum < 0.99)
        # find max gdist from effective sites
        # gdiff = gy - bkc.gcons[mi]
        # gmax = np.max(abs(gdiff))
        # distance.append(gmax)
        # find genetic distance enclosing middle 98% sites
        # dgmx, dgmn = gdiff.max(), gdiff.min()
        # gspan = dgmx - dgmn
        # distance.append(gspan)
        # get min, max positions in middle 98%
        # gmx, gmn = bkc.gcons[mi].max(), bkc.gcons[mi].min()
        # span = anno.gmap.interp_pos(gmx) - anno.gmap.interp_pos(gmn)
        # distance.append(span)
        # record the fraction of total sites contributing 98% of effect
        frac = 1.0 * np.sum(mi) / n
        fraction.append(frac)
        # record the conserved bases per Morgan in the middle 98%
        # b_g = np.sum(mi) / span
        # density.append(b_g)

    # plt.figure(figsize=(12, 7))
    # title = r'$\mathrm{proportion\ of\ conserved\ sites\ to\ reach\ 98\%\ ' \
    #         r'of\ b\ sum\ at\ a\ random\ sample\ of\ points}$'
    # plt.title(title, fontsize=18)
    # plt.xlabel(r'$\mathrm{proportion\ of\ sites}$', fontsize=18)
    # plt.xticks(fontsize=16)
    # plt.ylabel(r'$\mathrm{counts}$', fontsize=18)
    # plt.yticks(fontsize=16)
    # plt.hist(fraction, bins=100)

    # plt.figure(figsize=(15, 9))
    # plt.hist(distance, bins=100, label='max gdist')
    # plt.legend()
    #
    # plt.figure(figsize=(15, 9))
    # plt.hist(density, bins=100, label='density')
    # plt.legend()

    # plt.show()


def bkgd_range(pnt, cst, anno, coef, err):
    """find bounds for pos where remaining b terms sum to less than err"""

    # initialize BkgdCalculator for exact b calculations
    assert cst.chrom == anno.chrom
    udel = cst.fixed.u_fix
    bkc = BkgdCalculator(coef, udel, anno)

    # group sites into bins of uniform genetic distance
    g_start, g_end = bkc.gmap.gmap_limits
    g_bins = np.linspace(g_start, g_end, num=1e4)
    g_cnts = np.histogram(bkc.gcons, g_bins)[0]
    g_pnt = anno.gmap.interp_gpos(pnt)  # convert point to gmap units

    # find the index of the bin whose edges enclose pnt
    gi = np.searchsorted(g_bins, g_pnt)
    # calculate min genetic distance from g_pnt to each bin
    if gi < 2:
        # case 1: g_pnt either to the left of all bins or in the first bin
        g_dist = np.concatenate(([0], g_bins[1:-1] - g_pnt))
    elif gi > len(g_bins)-2:
        # case 2: g_ont is eiter right of all bins or in the last bin
        g_dist = np.concatenate((g_pnt - g_bins[1:-1], [0]))
    else:
        # case 3: g_pnt is in the middle with some bins on either side
        g_tuple = (g_pnt - g_bins[1:gi], [0], g_bins[gi + 1:] - g_pnt)
        g_dist = np.concatenate(g_tuple)

    # get bkgd effects on g_pnt from 1 single selected site at each bin dist
    b_bin = bkc.bkgd_sum(0.0, g_dist, bsum=False)

    # take the upper bound on the effects of bkgd from each bin by using
    # the minimum genetic distance to g_pnt for each binned site. this can
    # be done efficiently by multiplying the effects for a single site by the
    # by the conserved site count for the bin
    b_bin *= g_cnts

    # calculate the cumulative bkgd contribution per bin
    b_frac = np.cumsum(b_bin) / np.sum(b_bin)

    # use the cumulative contribution across bins to conservatively estimate
    # the bins indices accounting for >= (1-epsilon) of bkgd effects
    mi = np.searchsorted(b_frac, 0.5*err)-1
    mj = np.searchsorted(b_frac, 1-0.5*err)+1

    # cumulative conserved site COUNTS across bins become site INDICES, which
    # select the minimum set of conserved sites needed to make >= 1-epsilon
    # total b at g_pnt
    si = np.sum(g_cnts[:mi])
    sj = bkc.anno.num_bases - np.sum(g_cnts[mj:])

    return si, sj


def main():

    ch = 'chr20'
    anno = 'primate_cons95_Segments'
    coef = 10**-4.5
    err = 0.02

    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files, 0.01)
    ano = AnnoSegments(gmp, cst.bs_target(anno))

    # grid = bkgd_grid(cst, ano, cf, ep)
    # bterm_contribution(cst, ano, cf)
    # stop_sum(cst, ano, coef, err)


if __name__ == '__main__':
    main()
