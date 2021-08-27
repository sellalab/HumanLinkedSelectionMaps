import re
import sys
import numpy as np
from scipy.special import psi
from datetime import datetime as dt
from classes.runstruct import ChromStruct, root_dir, os
from classes.geneticmap import GeneticMap
from classes.annopoints import AnnoPoints
from classes.cscalculator import CSCalculator


__author__ = "davidmurphy"

"""
From the original MATLAB version of the program LS_PrecalcGridElements:
1. Get the neutral positions on the chromosome
2. Convert these to genetic map positions "gpos"
3. Thin the grid of neutral positions to look at based on the expected range at
which CS effects don't change for a given selection coefficient s
4. With this set of neutral positions in hand, estimate the coalescent rate due
to CS at substitutions for different values of s
"""

# effective pop sizes
ne_dmel = 1140800.0
ne_YRI = 2e4

# chromosome lengths dict
chr_length = dict(chr1=249250621, chr2=243199373, chr3=198022430,
                  chr4=191154276, chr5=180915260, chr6=171115067,
                  chr7=159138663, chr8=146364022, chr9=141213431,
                  chr10=135534747, chr11=135006516, chr12=133851895,
                  chr13=115169878, chr14=107349540, chr15=102531392,
                  chr16=90354753, chr17=81195210, chr18=78077248,
                  chr19=59128983, chr20=63025520, chr21=48129895,
                  chr22=51304566, chrY=59373566, chrX=155270560,
                  chr2L=23011544, chr2R=21146708, chr3L=24543557,
                  chr3R=27905053)


def ineq_rule(r1, r2, cmp_type=1):
    """
    Function used to determine if two coalescent rates can be distinguished
    from one another
    :param r1: the first rate
    :param r2: the second rate
    :param cmp_type: the type of comparison to use to perform the cutoff
    -- 1: the .3e rates from the standard .sw files are compared
    -- 2: round to .2e rates and compare
    -- 3: use .3e log(rates) to compare
    -- 4: use .2e log(rates) to compare
    -- 5: check if the ratio of the two rates is between 0.99 and 1.01
    :return: True if the two rates are judged to be different from one another,
    False if they are the same
    """
    if cmp_type is 1:
        return '{:.3e}'.format(r1) != '{:.3e}'.format(r2)
    elif cmp_type is 2:
        return '{:.2e}'.format(r1) != '{:.2e}'.format(r2)
    elif cmp_type is 3:
        r1, r2 = max(r1, 1e-6), max(r2, 1e-6)
        return '{:.3e}'.format(np.log(r1)) != '{:.3e}'.format(np.log(r2))
    elif cmp_type is 4:
        r1, r2 = max(r1, 1e-6), max(r2, 1e-6)
        return '{:.2e}'.format(np.log(r1)) != '{:.2e}'.format(np.log(r2))
    elif cmp_type is 5:
        r1, r2 = max(r1, 1e-6), max(r2, 1e-6)
        return not 0.99 < float(r1 / r2) < 1.01


def remove_duplicates(positions, rates, cmp_type=1):
    """
    STEP 1: remove consecutive duplicate rates -- if several lines have the
    same rate, keep only the HIGHEST position
    with the same rate
    :param positions: physical positions of the .sw map
    :param rates: coalescent rates for the .sw map
    :param cmp_type: the type of comparison to perform to determine duplicate
    rates (see function ineq_rule)
    :return pos, r: numpy arrays of the unique (consecutive) coalescent rates
    and highest associated position with each
    """
    pos, r = [positions[0]], [rates[0]]
    for i in xrange(1, len(rates)):
        # only keep consecutive sites when the rate changes within the rules
        if ineq_rule(r1=rates[i], r2=r[-1], cmp_type=cmp_type):
            r.append(rates[i])
            pos.append(positions[i])
        else:
            continue
    # extend the last pos to chr_length if it was cut off as a duplicate
    if pos[-1] < positions[-1]:
        pos[-1] = positions[-1]

    return map(np.array, [pos, r])


def interp_segments(positions, rates):
    """
    # STEP 2: find max segments where interpolation yields the same value
    (within error range)
    :param positions: physical positions (thinned of duplicates)
    :param rates: coalescent rates (thinned of duplicates)
    :return seg, r: numpy arrays of minimum segments and rates of coalescence
    in those segments
    """
    # convert positions to segments between each point
    segments = positions - np.concatenate([[0], positions[:-1]])
    seg, r = [segments[0]], [rates[0]]
    for i in xrange(1, len(rates)):
        # check 2 sample points between rate_1 and rate_2: 25% and 75%
        # from rate_1 to rate_2
        rate_1, rate_2 = r[-1], rates[i]
        delta_r = rate_2 - rate_1
        delta_x = segments[i]
        delta_rx = delta_r / delta_x
        # see if there is a detectable difference in the rate between
        # rate_1 and 50% to rate_2, starting from 1%
        fraction = 0.01
        size = 0
        while fraction < 0.5:
            delta_rate_i = rate_1 + fraction * delta_r
            if ineq_rule(rate_1, delta_rate_i, cmp_type=5):
                size = fraction * delta_x
                break
            else:
                fraction *= 2
        if size:
            steps = int(delta_x / size)
            size = int(size)
        else:
            steps = 0
        # for each 1% increase a new segment is created
        if steps > 1:
            for n in xrange(1, steps):
                # reduce the full segment by each inner segment
                segments[i] -= size
                seg.append(size)
                # interpolate the rate for inner segments
                interp_r = rate_1 + n * size * delta_rx
                r.append(interp_r)
            if segments[i]:
                seg.append(segments[i])
                r.append(rates[i])
        else:
            seg.append(segments[i])
            r.append(rates[i])

    return np.column_stack([r, seg])


def compress_csmap(cmap_file, cmp_type=1):
    """
    Take a pre-calculated CS map and reduce it further where possible
    (where positions share the same coalescent rate values within some level
    of error tolerance). Re-format the file to McVicker style segment lengths
    in the coordinate column for easier combination with the existing bmaps
    :param cmap_file: .npy file from generate_sw_maps
    :param cmp_type: the type of inequality rule to apply to rates to determine
    if they are the same within our
    precision levels
    :return compressed_map: compressed CS map
    """
    # get chrom for chr_length
    chrom = re.search('chr[\d\w]{1,2}', cmap_file).group()
    # load positions and rates into separate vectors
    if cmap_file.endswith('.npy'):
        pos, r = np.load(cmap_file).T
    else:
        pos, r = np.loadtxt(cmap_file).T
    assert pos[-1] <= chr_length[chrom]
    # fill out the last segment to extend to the end of the chrom
    pos[-1] = chr_length[chrom]
    # STEP 1: remove consecutive duplicate rates
    pos, r = remove_duplicates(positions=pos, rates=r, cmp_type=cmp_type)
    # STEP 2: find max segments where interpolation yields the same value
    # (within error range)
    # NOTE: this step makes the files 10x bigger! Does not work for compression!
    compressed = interp_segments(positions=pos, rates=r)

    return compressed


def cmap_header(cst, anno, coef, n_e):
    """Generate a header string for a sweep map file using ChromStruct"""
    sstr = '{:.8f}'.format(coef)
    subs_file = '{root}/data/csanno/{a}/substitutions_{chrom}_{a}.txt'
    gmap_file = '{root}/data/maps/{gmap}/genmap_{gmap}_{chrom}.txt'
    save_name = '{root}/precalc/{cmap_dir}/{gmap}_{chrom}_{a}_s{s}.sw'
    h = ['#{} {}'.format(cst.chrom, cst.chlen),
         '#CHROMOSOME_NAME={}'.format(cst.chrom),
         '#CHROMOSOME_LENGTH={}'.format(cst.chlen),
         '#OUTPUT_FILE=' + save_name.format(a=anno, s=sstr, **cst.dict),
         '#OUTPUT_TOKEN={}'.format(anno),
         '#RECOMB_RATE_TABLE=' + gmap_file.format(**cst.dict),
         '#RECOMB_RATE_SCALE=1.000000e-08',
         '#FOCALS_TABLE=' + subs_file.format(a=anno, **cst.dict),
         '#PARAM_NE0={:.6e}'.format(n_e),
         '#PARAM_S={s}'.format(s=sstr),
         '#PARAM_S_DIST_TYPE=POINT',
         '#TTF_APPROXIMATION=diffusion',
         '#StopSum=1',
         '#INTERP_METHOD=linear',
         '#MAX_DIST_SCALED=1',
         '#MAX_DIST=1.000000e+00']

    return '\n'.join(h)


def main_local():
    ch = 'chr20'
    cst = ChromStruct(ch, focal='nonsyn')
    fixed_precision_swmap(ch, cst.focal_files, cst.gmap_files, 0.01,
                          ne_YRI, 0.01)

    # generate_sw_maps(cst.chrom, cst.neut_masks, cst.focal_files,
    #                  cst.gmap_files,
    #                  coef=1e-3, ne=ne_YRI)
    # svals = 10**np.arange(-1.5, -6, -1)
    # fmt = '%d\t%.3e'
    # root = root_dir
    # fsubs = '{r}/data/csanno/{a}/substitutions_{c}_{a}.txt'
    # fsave = '{r}/precalc/dmel/comeron_{c}_{a}_s{s:.8f}.sw'
    # for ch in '2R 3L 3R'.split():
    #     neut = '{}/data/neut/dmel/Filtered_140521_poly_{}.txt'
    # .format(root, ch)
    #     gmap = '{}/data/maps/comeron/genmap_comeron_{}.txt'.format(root, ch)
    #     cst = ChromStruct(chrom=ch, gmap='comeron', cmap_dir='dmel')
    #     for an in 'UTR exonicNS intronic'.split():
    #         subs = fsubs.format(r=root, a=an, c=ch)
    #         for sx in svals:
    #             save = fsave.format(r=root, c=ch, a=an, s=sx)
    #             head = cmap_header(cst, an, sx, ne_dmel)
    #             sw = generate_sw_maps(neut, subs, gmap, sx, ne_dmel)
    #             np.savetxt(save, sw, fmt=fmt, header=head, comments='')


def main_remote_dmel():
    if len(sys.argv) != 4:
        print 'usage: sweep_maps chrom anno sval'
        exit(1)
    # record function start time
    start_time = dt.now()

    # read command line args
    ch, an = sys.argv[1:3]
    sx = eval(sys.argv[3])
    cst = ChromStruct(chrom=ch, gmap='comeron', cmap_dir='dmel')
    head = cmap_header(cst, an, sx, ne_dmel)
    fmt = '%d\t%.3e'
    root = root_dir

    # format input/output paths
    fneut = '{}/data/neut/dmel/Filtered_140521_poly_{}.txt'.format(root, ch)
    fgmap = '{}/data/maps/comeron/genmap_comeron_{}.txt'.format(root, ch)
    fsubs = '{r}/data/csanno/{a}/substitutions_{c}_{a}.txt'.format(
            r=root, a=an, c=ch)
    fsave = '{r}/precalc/dmel/comeron_{c}_{a}_s{s:.8f}.sw'.format(
            r=root, c=ch, a=an, s=sx)

    # get sweep map and save
    sw = generate_sw_maps(fneut, fsubs, fgmap, sx, ne_dmel)
    np.savetxt(fsave, sw, fmt=fmt, header=head, comments='')

    # print map name and total computation time
    etime = dt.now() - start_time
    print 'elapsed_time={}'.format(etime)


def diffusion_fixation_time(s, ne):
    """
    the diffusion approximation for fixation time (tau) of a beneficial allele
    starting from a new mutation with selection coefficient 's' and effective
    population size 'ne'
    :param s: the selection coefficient
    :param ne: effective population size
    :return tau: the diffusion approximation for tau:
    """
    euler = -psi(1)  # euler constant gamma
    return 2 * (np.log(4 * ne * s) + euler - 1 / (4 * ne * s)) / s


def probability_to_coalesce(genetic_distance, tau):
    """
    the probability that a site fixes due to a sweep at a nearby substitution
    :param genetic_distance: a vector of genetic distances from a putative
    sweep site
    :param tau:
    :return prob:
    """
    return np.exp(-genetic_distance * tau)


def load_data(gmp, neut_file, subs_file):
    """load neutral and focal pos, get gpos with GeneticMap instance"""
    # load all neutral sites from zipped 012 mask array
    assert neut_file.endswith('.npz')
    neut_pos, = np.where(np.load(neut_file)['neutmask'] != 0)

    # mask neutral sites beyond genetic mask and get their gpos
    neut_pos = gmp.gmap_mask(neut_pos)
    neut_gpos = gmp.interp_gpos(neut_pos)

    # gmap_mask substitution positions, convert to genetic map positions
    if subs_file.endswith('.npy'):
        sub_pos = np.load(subs_file)
    else:
        sub_pos = np.loadtxt(subs_file)

    # mask substitution sites beyond genetic mask and get their gpos
    sub_pos = gmp.gmap_mask(sub_pos)
    sub_gpos = gmp.interp_gpos(sub_pos)

    return neut_pos, neut_gpos, sub_gpos


def thin_grid(gdelta, gpos, glim):
    """
    Create a 'grid' of points between which the effects of linked selection
    with a given selection coefficient are indistinguishable. Find the nearest
    neutral sites to these points (if any) that can then be used to model the
    effects of sweeps given a set of substitutions. The effects at all other
    sites can be interpolated based on the predictions made at the thinned set
    of positions.
    :param gdelta: gdelta = the min distance where we might expect sweep
    effects to differ. below this distance (in genetic map units), we expect
    neutral sites to be effected by a sweep in the same way, so we can just
    focus on a single neutral site in these windows
    :param gpos: position of every neutral site in Morgans
    :param glim: the lower and upper bounds of the genetic map in Morgans
    :return grid_idx: indices for neutral sites where the effects of sweeps
    WILL become distinguishable given the
    selection coefficient and genetic map
    """
    # the spacing of genetic map intervals between which the effects of sweeps
    # are expected to be uniform. we can ignore all neutral sites except the
    # two closest to each spacing element (upstream and downstream)
    gmin, gmax = glim
    spacing = np.arange(gmin - gdelta, gmax + gdelta, gdelta)

    # remove duplicate gpos, keep the index, which correspond to the original
    # full size vector of positions. NOTE: unique sorts input, the return
    # indices are the indices that will produce the unique sorted array from
    # the original input array
    unique_gpos, unique_idx = np.unique(gpos, return_index=True)
    neut_edges = np.concatenate([[0], unique_gpos, [gmax + 1]])

    # we use hist to thin the grid:
    # some neutral sites are very densely packed and the effect of a
    # substitution will be the same for sites that are very close to one
    # another. thus we take a "grid" based on the selection coefficient that
    # is a rough estimate of the spacing at which the effects on neutral sites
    # will be indistinguishable. this allows us to save a much smaller
    # "pre-calculated" map of the effects of CS on select neutral sites, rather
    # than ALL neutral sites, and then later we can extrapolate the expected
    # effects on ALL sites during the inference from this condensed map
    neut_cnts, _ = np.histogram(a=spacing, bins=neut_edges)

    # find the bins where one or more of the spacing elements landed a hit
    # concatenate neutral position 'edges' for each pair of neutral sites that
    # capture one or more spacing elements hit_idx is scaled to unique_gpos,
    # unique_idx. make hit_idx unique and use to sort the original unique_idx
    # to get indices for neutral pos/gpos corresponding to the thinned points
    hits, = np.where(neut_cnts > 0)
    hit_idx = np.concatenate([hits[:-1], hits[1:]-1, [unique_gpos.size-1]])
    hit_idx = np.unique(hit_idx)
    grid_idx = unique_idx[hit_idx]

    # thin_gpos, thin_idx = unique_gpos[hit_idx], unique_idx[hit_idx]
    # thin_gpos, idx = np.unique(thin_gpos, return_index=True)
    # grid_idx = thin_idx[idx]

    return grid_idx


def fixation_rates(sub_gpos, tau, neut_gpos, s):
    """
    Iterate over the grid points and determine if any substitutions are close
    enough to affect fixation rates in the vicinity
    :param sub_gpos: genetic position of substitutions
    :param tau: diffusion approximate fixation time
    :param neut_gpos: neutral genetic positions
    :param s: selection coefficient
    :return: vector of fixation rate point estimates for each
    """
    # empty zeros array for fixation rates calculated at thinned neutral gpos
    fix_rates = np.zeros(neut_gpos.size)
    for (i, ngp) in enumerate(neut_gpos):
        # for each grid point look for any substitutions within "s" Morgans
        sidx,  = np.where(np.abs(sub_gpos - ngp) < s)
        if len(sidx):
            # calculate genetic distance to each sub (bounded by min_dist)
            dist = np.abs(sub_gpos[sidx] - ngp)
            # calculate probability to fix based on distance from substitution
            coal_prob = np.sum(probability_to_coalesce(dist, tau))
            fix_rates[i] = coal_prob

    return fix_rates
    # return np.column_stack([neut_pos[gidx], fix_rates])


def generate_sw_maps(chrom, neut_file, subs_file, gmap_file, coef, ne):
    """
    This function generates a map of the effects of selective sweeps (in terms
    of coalescent rate at a given site) for a set of neutral sites and a given
    selection coefficient. Map is saved to an output file directly
    :param chrom: chromosome name
    :param neut_file: neutral sites
    :param subs_file: substitutions (physical pos)
    :param gmap_file: genetic map
    :param coef: selection coefficient s
    :param ne: the effective population size for the current sample
    """
    # calculate tau once per coefficient
    tau = diffusion_fixation_time(coef, ne)

    # average genetic distance between two bases for map
    # equivalent to 3 bp using the Hinch et al., 2011 African American Map:
    # - genetic distance across autosomes = 35.234 Morgans
    # - bases across autosomes = 2881033286 bp
    # - minimum_neutral_distance = 3bp * 35.234M / 2881033286bp ~ 4e-08
    minimum_grid_distance = 3e-8
    grid_delta = max(coef * 0.01, minimum_grid_distance)

    # initialize GeneticMap instance
    gmp = GeneticMap(chrom, gmap_file, gscale=0.01)

    # load neutral pos/gpos and subs in gpos
    neut_pos, neut_gpos, sub_gpos = load_data(gmp, neut_file, subs_file)

    # use indices of maximum CS effect resolution for "coef" to thin neut sites
    gidx = thin_grid(grid_delta, neut_gpos, gmp.gmap_limits)
    neut_pos = neut_pos[gidx]
    neut_gpos = neut_gpos[gidx]

    # create array of [neut_pos, fixation_rate] for each grid point
    fix_rates = fixation_rates(sub_gpos, tau, neut_gpos, coef)

    return np.column_stack((neut_pos, fix_rates))


def fixed_precision_swmap(ch, fsubs, fgmap, coef, ne, err, lim=True):
    """calculate the effects of sweeps at each point on a fixed grid"""
    # minimum grid interval is approximately 1bp in Morgans
    min_grid_interval = 1.22e-08

    # load data classes
    gmp = GeneticMap(ch, fgmap, 0.01)
    ano = AnnoPoints(gmp, afile=fsubs)
    csc = CSCalculator(ano, coef, min_grid_interval, ne)

    # use grid in Morgans for CS calculations, use bp grid for final table
    grid = csc.grid_gpos(err)
    grid_pos = gmp.interp_pos(grid)

    # use CSCalculator to calculate the effects of CS at each  grid point
    fix_rates = csc.grid_fixrates(lim=lim)
    # fix_rates = csc.fixrate_vector(grid, gpos=True, lim=lim)

    return np.column_stack((grid_pos, fix_rates))


def main_remote_human():
    if len(sys.argv) != 5:
        print 'usage: sweep_maps chrom anno sval lim'
        exit(1)

    # record function start time
    start_time = dt.now()

    # fixed params
    err = 0.01
    ne = 2e4  # Ne for YRI
    tk = 'pr95.cleanrun'  # standard clean run data token

    # read command line args
    ch = sys.argv[1]
    an = sys.argv[2]
    coef = eval(sys.argv[3])
    lim = eval(sys.argv[4])

    # load chrom struct
    cdir = '{}_fixed_cmap'.format(an)
    cst = ChromStruct(ch, cmap_dir=cdir, tkn=tk, focal=an)
    # create the new map folder if it does not yet exist
    cpath = '{root}/precalc/{cmap_dir}'.format(**cst.dict)
    if not os.path.isdir(cpath):
        os.mkdir(cpath)

    # file save paths
    fsubs = cst.focal_files
    fgmap = cst.gmap_files
    fsave = cst.cmap_file(an, coef)
    # edit file suffix if no limit is used
    if not lim:
        fsave = fsave.replace('{}.s0'.format(an), '{}.nolim.s0'.format(an))

    # get new fixed grid map, save copy in numpy file format
    sw = fixed_precision_swmap(ch, fsubs, fgmap, coef, ne, err, lim)
    np.save(fsave.replace('.sw', '.npy'), sw)

    # sw = generate_sw_maps(ch, fneut, fsubs, fgmap, coef, ne)
    # # format the sweep map file for MATLAB readability
    # head = cmap_header(cst, an, coef, ne_dmel)
    # fmt = '%d\t%.3e'
    # np.savetxt(fsave, sw, fmt=fmt, header=head, comments='')

    # print map name and total computation time # get sweep map and save
    etime = dt.now() - start_time
    print 'elapsed_time={}'.format(etime)


if __name__ == "__main__":
    if root_dir.startswith('/Users/davidmurphy'):
        main_local()
        # main_remote_human()
        # pass
    else:
        main_remote_human()
