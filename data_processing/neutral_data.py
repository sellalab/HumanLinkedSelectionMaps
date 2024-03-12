import numpy as np
import data_tools as dtools
from functions import calc_pi
from classes.bkgdmap import BkgdMap
from sys import argv, stderr, stdout
from classes.runstruct import ChromStruct, default_filter_file, root_dir

xdir = root_dir.replace('linked_selection', 'x_aut')
fpth = 'NONE'

__author__ = 'davidmurphy'

# set token externally
# token = 'xaut_rmcpgi'
# token = 'x_aut_ratios'
# token = 'xaut_rmbgc'


class FilePaths:
    """structure organizing file paths and directories used for X/A analyses"""
    def __init__(self, token):
        # token for filepaths
        self.token = token

        # XA analysis path
        self.root = root_dir.replace('linked_selection', 'x_aut')

        # aut bval file path, X & A b-index file paths
        self.abval_file = '{}/joined/aut.bval.{}.npz'.format(self.root, token)
        self.anewbval_file = '{}/joined/aut.newbval.{}.npz'.format(self.root,
                                                                   token)
        self.xbidx_file = '{}/x.bidx.{}.npz'.format(self.root, token)
        self.bidx_file = '{}/joined/aut.bidx.{}.npz'.format(self.root, token)

        # aut dist file path, X & A dist index file paths
        self.adist_file = '{}/joined/aut.dist.{}.npz'.format(self.root, token)
        self.xdidx_file = '{}/x.didx.{}.npz'.format(self.root, token)
        self.didx_file = '{}/joined/aut.didx.{}.npz'.format(self.root, token)

    def fbval(self, chrom):
        """bval file paths for individual chrom"""
        f_str = '{}/bvals/{}.bvals.{}.npz'
        return f_str.format(self.root, chrom, self.token)

    def fdist(self, chrom):
        """dist file paths for individual chrom"""
        return '{}/dists/{}.mindists.{}.npz'.format(self.root, chrom,
                                                    self.token)

    def fpi(self, chrom, neut):
        """pi file paths for individual chrom and population"""
        return '{}/pivals/{}.{}.pi.{}.npz'.format(self.root, chrom, neut,
                                                  self.token)

    def fapi(self, neut):
        """joined aut pi values for specific population (=neut)"""
        return '{}/joined/aut.{}.pi.{}.npz'.format(self.root, neut, self.token)

    def fdiv(self, chrom, outg):
        """ div file paths for individual chrom and outgroup"""
        return '{}/dvals/{}.{}.d.{}.npz'.format(self.root, chrom, outg,
                                                self.token)

    def fadiv(self, outg):
        """joined aut div values for specific population (=neut)"""
        return '{}/joined/aut.{}.div.{}.npz'.format(self.root, outg, self.token)

    def fbinpi(self, neut, sortby):
        """binned X & A pi values"""
        return '{}/saves/{}.pi.s{}.{}.npz'.format(self.root, neut, sortby,
                                                  self.token)

    def fbindiv(self, outg, sortby):
        """binned X & A div values"""
        return '{}/saves/{}.div.s{}.{}.npz'.format(self.root, outg, sortby,
                                                   self.token)

    def fbincnt(self, sortby):
        """bins of sort values, X & A neutral site counts per bin"""
        return '{}/saves/cnts.s{}.{}.npz'.format(self.root, sortby, self.token)

    def fblck(self, neut, sortby, bnum, mcv=True):
        if mcv:
            fstr = '{}/blocks/{}.blck.s{}.n{}.{}.npz'
        else:
            fstr = '{}/blocks/{}.blck.snew{}.n{}.{}.npz'
        return fstr.format(self.root, neut, sortby, bnum, self.token)

    def fblck_list(self, neut, sortby, wsize):
        fstr = '{}/blocks/{}.blklist.s{}.wsize_{:.2e}.{}.npz'
        return fstr.format(self.root, neut, sortby, wsize, self.token)

    def fxblck_list(self, neut, sortby, wsize):
        fstr = '{}/blocks/{}.xblklist.s{}.wsize_{:.2e}.{}.npz'
        return fstr.format(self.root, neut, sortby, wsize, self.token)

    def fxblck(self, neut, sortby, bnum):
        fstr = '{}/blocks/{}.xblck.s{}.n{}.{}.npz'
        return fstr.format(self.root, neut, sortby, bnum, self.token)

    def fpibd(self, neut, bmin, bmax):
        fstr = '{}/resample_B/{}.pi_bd.b{}_{}.{}.npz'
        return fstr.format(self.root, neut, bmin, bmax, self.token)

    def fxpibd(self, neut, bmin, bmax):
        fstr = '{}/resample_B/{}.xpi_bd.b{}_{}.{}.npz'
        return fstr.format(self.root, neut, bmin, bmax, self.token)

    def fpid(self, neut, thresh, wsize):
        fstr = '{}/resample_cM/{}.pi_d.thresh.{}.cM.wsize_{:.2e}.{}.txt'
        return fstr.format(self.root, neut, thresh, wsize, self.token)

    def fxpid(self, neut, thresh, wsize):
        fstr = '{}/resample_cM/{}.xpi_d.thresh.{}.cM.wsize_{:.2e}.{}.txt'
        return fstr.format(self.root, neut, thresh, wsize, self.token)

    def f_allblk(self, neut, xaut, bmin, bmax):
        fstr = '{}/saves/{}.{}.allblocks.b{}_{}.{}.txt'
        return fstr.format(self.root, neut, xaut, bmin, bmax, self.token)

# load struct of file paths

# biased gene conversion pairs
is_bgc = ['AC', 'AG', 'CA', 'CT', 'GA', 'GT', 'TC', 'TG']
non_bgc = ['AA', 'AT', 'CC', 'CG', 'GC', 'GG', 'TA', 'TT']


def get_neutral_masks(cst):
    """
    A variation on neutrality_mask function in data_tools. Two neutral site
    masks (one encoding macaque divergence, the other encoding orangutan
    divergence) are constructed in parallel and as a final step sites missing
    data for either species are filtered from each mask.
    """
    # phastcons 100-vert p=0 sites
    consmask = dtools.conservation_mask(cst.wigz_ncons, cst.chrom)

    # strict calling mask from 1000 genomes project
    callmask = dtools.sequence_mask(cst.call_mask, cst.chrom)

    # aligned positions to macaque (with substitutions set to 2)
    if cst.tkn == 'xaut_rmbgc':
        submask = dtools.substitution_mask(cst.axt_outgroup, cst.chrom, is_bgc)
    elif cst.tkn == 'xaut_kpbgc':
        submask = dtools.substitution_mask(cst.axt_outgroup, cst.chrom, non_bgc)
    else:
        submask = dtools.substitution_mask(cst.axt_outgroup, cst.chrom)

    # aggregate filter regions into a binary mask array
    filtmask = dtools.get_filter_mask(cst.chrom, default_filter_file)

    # SNP sites (used to encode polymorphism in the final neutrality mask
    fsnp = cst.snp_files.replace('.phase3', '.female.phase3')
    # remove BGC or non-BGC via the filtermask at BGC SNPs
    if cst.tkn == 'xaut_rmbgc':
        snpcnt = dtools.snpcount(fsnp, returnbases=True)
        issnp, ref, alt = snpcnt[0], snpcnt[1], snpcnt[3]
        pairs = np.core.defchararray.add(ref, alt)
        ibgc = np.in1d(pairs, is_bgc)
        filtmask[(issnp-1)[ibgc]]= 0
    elif cst.tkn == 'xaut_kpbgc':
        snpcnt = dtools.snpcount(fsnp, returnbases=True)
        issnp, ref, alt = snpcnt[0], snpcnt[1], snpcnt[3]
        pairs = np.core.defchararray.add(ref, alt)
        ibgc = np.in1d(pairs, non_bgc)
        filtmask[(issnp - 1)[ibgc]] = 0
    else:
        issnp = dtools.snpcount(fsnp)[0]

    # add CpG sites to the filter mask (ref and complement strands)
    if cst.tkn != 'xaut_rmcpgi':
        anc_fname = cst.ancs_files
        cpg = dtools.cpg_mask(cst.chrom, anc_fname)
        filtmask[cpg] = False
        filtmask[cpg+1] = False

    # add CpG islands to filter mask
    cpgi = np.loadtxt(cst.cpg_isld, usecols=(1, 2), dtype=int)
    filtmask = dtools.mask_segments(filtmask, cpgi)

    # remove chrX PAR (source: http://genome.ucsc.edu/cgi-bin/hgGateway)
    if cst.chrom == 'chrX':
        filtmask[60000:2699520] = False
        filtmask[154931043:155260560] = False

    # neutral mask encoding macaque divergence
    nmsk_args = (consmask, callmask, submask, issnp, filtmask, None)
    nmsk = dtools.neutrality_mask(*nmsk_args)

    return nmsk


def nearest_feature(feature_map, point_map, genetic_map):
    """
    for each point in point_map, find the nearest feature in feature_map
    :param feature_map: a set of coordinates such as  start/end of exons
    :param point_map: a set of genomic positions
    :param genetic_map: genetic map used to convert features and points into
    genetic map units
    :return: feature_distances: distance to nearest feature for each point
    in genetic map units
    """

    # use the genetic map to convert features and points
    pos, gpos = genetic_map.T
    gi, gj = np.interp(feature_map, pos, gpos).T
    gpoints = np.interp(point_map, pos, gpos)

    # use the unique set of genetic feature minpoints to stratify the points
    gbins = np.unique(gi + 0.5*(gj-gi))

    # get the indices of the lower and upper bin edges for each point
    jdx = np.minimum(np.searchsorted(gbins, gpoints), len(gbins)-1)
    idx = np.maximum(0, jdx-1)

    # min distance to upper and lower bin edge = nearest feature distance
    dists = abs(np.column_stack((gbins[jdx] - gpoints, gpoints - gbins[idx])))
    min_dists = np.min(dists, axis=1)

    return min_dists


def mcvicker_b(point_map, bmap):
    # load bvals scaled from 1-1000 into BkgdMap
    bvals, segs = bmap.T
    bmp = BkgdMap(bvals, segs)

    # find bmap segments and corresponding b estimates for each neutral point
    b = bmp.interp_bkgd(point_map)

    return b


def pi_vals(cst, neut):
    """calculate pi for each polymorphic neutral site"""
    # load SNP count data
    fsnp = cst.snp_files.replace('.phase3', '.female.phase3')
    snpcnt = np.column_stack(dtools.snpcount(fsnp))

    # calculate pi for subset of neutral SNPs
    neutsnps = dtools.neutral_snps(neut, snpcnt)
    pos, sample, alt = neutsnps.T
    pi = calc_pi(sample+alt, alt)

    # reduce neut to the set of valid neutral sites
    neut = neut[neut != 0]

    # record pi for each neutral site (default = 0)
    neutpi = np.zeros(neut.size)
    neutpi[(neut > 2)] = pi

    return neutpi


def d_vals(neut):
    """get divergence status for each neutral site"""
    neut = neut[neut != 0]
    div = np.zeros(neut.size)
    div[((neut == 2) | (neut == 6))] = 1

    return div


def save_pivals(cst, nmsk, fpth):
    """save pi and div at each neutral site as compressed arrays"""
    # calculate pi across neutral sites
    pi = pi_vals(cst, nmsk)

    # saved compressed neutral pi values
    np.savez_compressed(fpth.fpi(cst.chrom, cst.neut), pi=pi)
    stderr.write('{} neutral pi values saved\n'.format(cst.chrom))
    stdout.flush()


def save_dvals(cst, nmsk, fpth):
    """save div at each neutral site (sub=1, match=0)"""
    # get array of divergence for neutral sites, with subs=1
    div = d_vals(nmsk)

    # saved compressed divergence values
    np.savez_compressed(fpth.fdiv(cst.chrom, cst.outg), div=div)

    msg = '{chrom} human:{outg} D values saved\n'.format(**cst.dict)
    stderr.write(msg)
    stdout.flush()


def save_mindists(cst, nmsk, fpth):
    """save minimum genetic distance to exon for each neutral site"""
    feats = np.loadtxt(cst.nr_exons, usecols=(1, 2))

    # load genetic map
    gmap = np.loadtxt(cst.gmap_files, usecols=(0, 2))

    # sex average chrX map
    if cst.chrom == 'chrX':
        gmap[:,1] *= 2.0/3.0

    # load neutral site positions
    pts = np.where(nmsk != 0)[0]

    # calculate the distance to the nearest feature for each neutral point
    dist = nearest_feature(feats, pts, gmap)

    # save compressed distance values for each neutral site
    np.savez_compressed(fpth.fdist(cst.chrom), dist=dist)

    stderr.write('{} cM to nearest exon saved\n'.format(cst.chrom))
    stdout.flush()


def save_bvals(cst, nmsk, fpth):
    """save McVicker B at each neutral site"""
    # load B maps
    bfile = '{root}/mcvicker/hg19_bkgd/{chrom}.hg19.bkgd'.format(**cst.dict)
    bmap = np.loadtxt(bfile)

    # load neutral site positions
    pts = np.where(nmsk != 0)[0]

    # interpolate b values for neutral sites
    bval = mcvicker_b(pts, bmap)

    # save compressed b vals for each neutral site
    np.savez_compressed(fpth.fbval(cst.chrom), bval=bval)

    stderr.write('{} neutral site B values saved\n'.format(cst.chrom))
    stdout.flush()


def get_arrays(chrom, neut, outg, token):
    """get pi, div, dist, bvals for each chrom one at a time"""
    # create the ChromStruct
    cst = ChromStruct(chrom=chrom, tkn=token, neut=neut, outg=outg)
    # initialize filepaths data struct
    fpth = FilePaths(token=cst.tkn)

    # generate the neutral mask used by functions
    nmsk = get_neutral_masks(cst)

    # save pi values for each different population
    save_pivals(cst, nmsk, fpth)

    # only save div, dists and bvals once -- these are the same for all pops
    if neut == 'YRI':
        save_dvals(cst, nmsk, fpth)
        save_mindists(cst, nmsk, fpth)
        save_bvals(cst, nmsk, fpth)

    return None


def join_arrays(neut, outg, token, x=0, pi=0, div=0, b=0, dist=0):
    """join pi, div, dists & bvals across autosomes"""
    cst = ChromStruct(chrom='chr1', tkn=token, neut=neut, outg=outg)
    # initialize filepaths data struct
    fpth = FilePaths(token=cst.tkn)

    """
    AUTOSOME DATA
    """
    # collect lists of each file type to process one at a time
    fp, fd, fb, fm = [], [], [], []
    for ch in cst.chroms:
        cst.chrom = ch
        fp.append(fpth.fpi(ch, neut))
        fd.append(fpth.fdiv(ch, cst.outg))
        fb.append(fpth.fbval(ch))
        fm.append(fpth.fdist(ch))

    # concatenate and save aut pi data, delete array from memory
    if pi:
        dpi = np.concatenate([np.load(f)['pi'] for f in fp])
        pifile = fpth.fapi(cst.neut)
        np.savez_compressed(pifile, pi=dpi)
        del pi
        stderr.write('neutral aut pi values concatenated and saved\n')
        stdout.flush()

    # concatenate and save aut div data, delete array from memory
    if div:
        ddiv = np.concatenate([np.load(f)['div'] for f in fd])
        divfile = fpth.fadiv(cst.outg)
        np.savez_compressed(divfile, div=ddiv)
        del div
        stderr.write('neutral aut div values concatenated and saved\n')
        stdout.flush()

    # concatenate and save aut bval data
    if b:
        bval = np.concatenate([np.load(f)['bval'] for f in fb])
        bvalfile = fpth.abval_file
        np.savez_compressed(bvalfile, bval=bval)
        stderr.write('neutral aut B values concatenated and saved\n')
        stdout.flush()

        # get sorting indices for concatenated bvals and save
        bidx = np.argsort(bval)
        del bval  # delete bvals from memory after saving and using to sort
        bidxfile = fpth.bidx_file
        np.savez_compressed(bidxfile, bidx=bidx)
        del bidx  # delete bidx after saving
        stderr.write('neutral aut B value sorting index saved\n')
        stdout.flush()

    # concatenate and save aut exon distance data
    if dist:
        ddist = np.concatenate([np.load(f)['dist'] for f in fm])
        distfile = fpth.adist_file
        np.savez_compressed(distfile, dist=ddist)
        stderr.write('neutral aut exon distances concatenated and saved\n')
        stdout.flush()

        # get exon distance sorting indices for concatenated distances and save
        didx = np.argsort(dist)
        del ddist  # delete dist from memory after using to sort
        didxfile = fpth.didx_file
        np.savez_compressed(didxfile, didx=didx)
        del didx  # delete didx after saving
        stderr.write('neutral aut exon distances sorting index saved\n')
        stdout.flush()

    """
    X DATA
    """
    if x:
        # get sorting indices for chrX bvals
        fb_x = fpth.fbval('chrX')
        xbidx = np.argsort(np.load(fb_x)['bval'])
        bidxfile = fpth.xbidx_file
        np.savez_compressed(bidxfile, xbidx=xbidx)
        stderr.write('neutral chrX B value sorting index saved\n')
        stdout.flush()

        # get sorting indices for chrX exon distances
        fm_x = fpth.fdist('chrX')
        xdidx = np.argsort(np.load(fm_x)['dist'])
        didxfile = fpth.xdidx_file
        np.savez_compressed(didxfile, xdidx=xdidx)
        stderr.write('neutral chrX exon distances sorting index saved\n')
        stdout.flush()

    return None


def sort_idx(sortby='bval'):
    """get index to sort whole X & A neutral site arrays"""
    if sortby == 'bval':
        xidx = np.load(fpth.xbidx_file)['xbidx']
        aidx = np.load(fpth.bidx_file)['bidx']
    elif sortby == 'dist':
        xidx = np.load(fpth.xdidx_file)['xdidx']
        aidx = np.load(fpth.didx_file)['didx']
    else:
        raise KeyError('no files matching {}'.format(sortby))

    # indices within uint32 range so we can cut the memory use in half
    return xidx.astype('uint32'), aidx.astype('uint32')


def bin_index(nbins, sortby):
    """
    Get coordinates to split data into n bins of cM to exon or B value.

    Generate n uniformly spaced bins spanning the range of the feature to be
    binned. For example, if B values range from 1 to 1000 and nbins=100, define
    the bins as: [(1,10), (11, 20)...(981,990), (991, 1000)]).

    :param nbins: number of bins to sort data into
    :param sortby: statistic used for sorting -- cM from exon or B val
    :return:
    """
    # get sorting indices for x, aut files
    xidx, aidx = sort_idx(sortby=sortby)

    # load X/A bval/dist data sorted using b sorting index
    if sortby == 'bval':
        # load bvals (discrete valued 1-1000 so use uint16 to save memory)
        xval = np.load(fpth.fbval('chrX'))['bval'][xidx].astype('uint16')
        aval = np.load(fpth.abval_file)['bval'][aidx].astype('uint16')

        # bins are always 0-1000
        bin_vals = np.arange(1002, dtype='uint16')

    elif sortby == 'dist':
        # load distances
        xval = np.load(fpth.fdist('chrX'))['dist'][xidx]
        aval = np.load(fpth.adist_file)['dist'][aidx]

        # get max values for dist bins (min is set to 0)
        max_val = max(xval.max(), aval.max())

        # # create nbins uniformly spaced bins from min to max
        # bin_vals = np.linspace(start=min_val, stop=max_val, num=nbins)

        # use the bin size from Hernandez et al., 2011
        hndz_bin = 1.2e-05

        # bin up to the limit used in Hernandez
        # bin_vals = np.arange(0, 0.4+hndz_bin, hndz_bin)

        # bin through the maximum cM in the data
        # bin_vals = np.arange(0, max_val+hndz_bin, hndz_bin)
        cm_close = np.arange(0, 0.25125, 0.00125)
        cm_far = np.arange(0.25125, 1.01, 0.01)
        bin_vals = np.concatenate((cm_close, cm_far, [10]))

    else:
        # raise an error if the key is not recognized
        raise KeyError('no files matching {}'.format(sortby))

    # use searchsorted for the complete set of (start, end) indices
    bin_xidx = np.searchsorted(xval, bin_vals).astype('uint32')
    bin_aidx = np.searchsorted(aval, bin_vals).astype('uint32')

    # NOTE: cut first bin edge (0) so bin_vals length matches means and counts
    return xidx, aidx, bin_xidx, bin_aidx, bin_vals[1:]


def bin_means(bin_idx, xval, aval):
    """get mean X/A data values in bins delineated by idx"""
    # initialize empty arrays for X/A averages
    nbins = len(bin_idx)-1
    xavg, aavg = np.zeros(shape=nbins), np.zeros(shape=nbins)

    # get mean X/A values per bin
    for i in range(nbins):
        # get (start, end) indices for X/A
        x_start, a_start = bin_idx[i]
        x_end, a_end = bin_idx[i+1]

        # take the mean from start:end
        xavg[i] = np.mean(xval[x_start:x_end])
        aavg[i] = np.mean(aval[a_start:a_end])

    return xavg, aavg


def bin_pi(bin_idx, xidx, aidx, neut):
    """average sorted pi values in bins defined by idx"""
    # sort pi arrays as they're loaded
    xpi = np.load(fpth.fpi('chrX', neut))['pi'][xidx]
    api = np.load(fpth.fapi(neut))['pi'][aidx]

    # get binned averages for pi values
    avg_xpi, avg_api = bin_means(bin_idx, xval=xpi, aval=api)

    return avg_xpi, avg_api


def bin_div(bin_idx, xidx, aidx, outg):
    """average sorted D values in bins defined by idx"""
    # load sorted X/A div data
    xdiv = np.load(fpth.fdiv('chrX', outg))['div'][xidx]
    adiv = np.load(fpth.fadiv(outg))['div'][aidx]

    # get binned averages for div values
    avg_xdiv, avg_adiv = bin_means(bin_idx, xval=xdiv, aval=adiv)

    return avg_xdiv, avg_adiv


def bin_arrays():
    if root_dir.startswith('/Users/davidmurphy/'):
        nbins = 1000
        neut = 'YRI'
        outg = 'rheMac3'
        sortby = 'bval'

    else:
        nbins = int(eval(argv[1]))
        neut = argv[2]
        outg = argv[3]
        sortby = argv[4]

    # get sorting indices for X/A neutral data, binning indices and bin values
    xidx, aidx, xbidx, abidx, bins = bin_index(nbins, sortby)
    # combine binning indices into one array
    bidx = np.column_stack((xbidx, abidx)).astype('uint32')
    stderr.write('loaded X/A sorting indices and bin edges\n')
    stdout.flush()

    # get counts for each bin
    acnt = abidx[1:] - abidx[:-1]
    xcnt = xbidx[1:] - xbidx[:-1]

    # save counts and bins
    f = fpth.fbincnt(sortby)
    np.savez_compressed(f, bins=bins, acnt=acnt, xcnt=xcnt)
    stderr.write('saved X/A {} bins and site counts\n'.format(sortby))
    stdout.flush()

    # get X/A mean binned pi for each pop
    for neut in 'YRI LWK GIH CEU CHB TSI JPT MXL'.split():
        avg_xpi, avg_api = bin_pi(bidx, xidx, aidx, neut)
        f = fpth.fbinpi(neut, sortby)
        np.savez_compressed(f, api=avg_api, xpi=avg_xpi)
        stderr.write('saved mean X/A {} pi in bins\n'.format(neut))
        stdout.flush()

    # get X/A mean binned div values for each outgroup
    for outg in ['rheMac3']:
        avg_xdiv, avg_adiv = bin_div(bidx, xidx, aidx, outg)
        f = fpth.fbindiv(outg, sortby)
        np.savez_compressed(f, adiv=avg_adiv, xdiv=avg_xdiv)
        stderr.write('saved mean X/A {} div in bins\n'.format(outg))
        stdout.flush()

    return None


def gc_content(chrom, outg, bminmax=None, cmin=None):
    """calculate the GC content for neutral data with optional B/cM ranges"""
    # generate ChromStruct
    cst = ChromStruct(chrom=chrom, outg=outg)

    # create neutral mask and get neutral mask positions
    nmsk = get_neutral_masks(cst)
    nidx = np.where(nmsk > 0)[0]

    # if there are b or cm filters, apply to nidx
    if bminmax is not None:
        filt = '{}<=B<={}'.format(*bminmax)
        # bminmax is a tuple of minimum and maximum values to keep
        mn, mx = bminmax
        # generate the bvalues (don't used save as it may be rheMac or ponAbe)
        bfile = '{root}/mcvicker/hg19_bkgd/{chrom}.hg19.bkgd'.format(**cst.dict)
        bmap = np.loadtxt(bfile)
        bval = mcvicker_b(nidx, bmap)
        # create a mask for bvals in the minmax range
        bmsk = (bval >= mn) & (bval <= mx)
        nidx = nidx[bmsk]

    elif cmin is not None:
        filt = 'cM>={}'.format(cmin)
        # generate the dists (don't used save as it may be rheMac or ponAbe)
        feats = np.loadtxt(cst.nr_exons, usecols=(1, 2))
        gmap = np.loadtxt(cst.gmap_files, usecols=(0, 2))
        # sex average chrX map
        if cst.chrom == 'chrX':
            gmap[:, 1] *= 2.0 / 3.0
        # calculate the distance to the nearest feature for each neutral point
        dist = nearest_feature(feats, nidx, gmap)
        # create a mask for cm values above the cmin range
        cmsk = dist >= cmin
        nidx = nidx[cmsk]

    else:
        filt = 'NONE'

    # get ref genome sequence
    ref = dtools.fasta_array(chrom, cst.refg_files)

    # count ATs and CGs
    at = np.sum(np.in1d(ref[nidx], ['A', 'T']))
    cg = np.sum(np.in1d(ref[nidx], ['C', 'G']))

    # print info for the chromosome
    print('{} {} {} FILTER={}'.format(chrom, at, cg, filt))


def main():
    if root_dir.startswith('/Users/davidmurphy/'):

        ch, nt, og = ['chr22', 'TSI', 'rheMac3']
        bmnmx = (600, 900)
        cmin = 0.2
        chroms = ['chr{}'.format(c) for c in range(1,23)+['X']]
        # print 'B value filters:'
        # for c in chroms:
        #     gc_content(c, og, bminmax=bmnmx)
        print('\ncM value filters:')
        for c in chroms:
            gc_content(c, og, cmin=cmin)

    else:
        if len(argv) == 4:
            ch, og, tp = argv[1:4]
            if tp == 'bval':
                bmnmx = (600, 900)
                gc_content(ch, og, bminmax=bmnmx)
            elif tp == 'cmin':
                cmin = 0.2
                gc_content(ch, og, cmin=cmin)
            else:
                gc_content(ch, og)

        if len(argv) == 5:
            ch, nt, og, tk = argv[1:5]
            get_arrays(ch, nt, og, tk)
        if len(argv) == 9:
            nt, og, tk = argv[1:4]
            x, pi, div, b, dist = [eval(arg) for arg in argv[4:]]
            join_arrays(nt, og, tk, x, pi, div, b, dist)

if __name__ == '__main__':
    main()
