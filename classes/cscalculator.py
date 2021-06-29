import numpy as np
from itertools import izip
from scipy.special import psi
from annopoints import AnnoPoints
from geneticmap import GeneticMap
from datetime import datetime as dt
from sys import argv, stderr, stdout
from runstruct import ChromStruct, root_dir
from data_processing.data_tools import randint_unique, calculate_error


__author__ = 'davidmurphy'


class SweepMap(object):
    """a class to represent sweep maps and perform basic functions"""
    def __init__(self, cfile):
        self.cfile = cfile
        self.pos, self.css = self._readmap()

    def __len__(self):
        """defined as the number of segments in the B map"""
        return len(self.css)

    def _readmap(self):
        """the reader handles 3 different cs map formats"""
        if self.cfile.endswith('.npy'):
            return np.load(self.cfile).T
        elif self.cfile.endswith('.npz'):
            return np.load(self.cfile)['sweepmap'].T
        else:
            return np.loadtxt(self.cfile).T

    def interp_idx(self, pos):
        """get the sweep map index closest to position"""
        imax = len(self) - 1
        idx = np.minimum(np.searchsorted(self.pos, pos), imax)
        return idx

    def interp_cs(self, pos):
        """get cs effect at a given position on the map"""
        idx = self.interp_idx(pos)
        return self.css[idx]


class CSCalculator(object):
    """a class that calculates the effects of CS on neutral fixation rates"""
    def __init__(self, anno, coef, ne, err=0.01, rcalc='haldane'):
        assert isinstance(anno, AnnoPoints)
        self.anno = anno  # AnnoPoints instance containing substitutions
        self.coef = coef  # selection coefficient
        self.ne = ne  # effective population size
        self.err = err # relative error to enforce block boundaries
        self.rcalc = rcalc  # method for calculating rec between points

    def _rxy(self, gy, gx):
        """calculate linear distance or Haldane's r between gmap points"""
        m = abs(gy - gx)
        if self.rcalc == 'linear':
            return m
        elif self.rcalc == 'haldane':
            return -np.expm1(-2.0*m) / 2.0  # NOTE: expm1 protects precision
        else:
            raise ValueError(self.rcalc + ' not an acceptable rcalc method')

    def fixrate_point(self, gy, lim=False, rmin=False):
        """summed effect of CS on neutral fixation rate for site y"""
        # if gy is array or list, call recursively and return array
        if isinstance(gy, np.ndarray or isinstance(gy, list)):
            return np.array([self.fixrate_point(y) for y in gy])

        # make sure gy is a genetic map position inside of the currant gmap
        if not (self.gmap.gpos.min() <= gy <= self.gmap.gpos.max()):
            msg = '{} genetic map position {} not found'
            raise ValueError(msg.format(self.chrom, gy))
        # get rec probabilities between neutral point and sweep sites
        gdist = self._rxy(gy, self.anno.rpoints)
        # optional mask of subs >= s Morgans away from neutral site
        if lim:
            gdist = gdist[(gdist < self.coef)]

        # get coalescence probabilities at y accounting for sweeps across anno
        probs = self.coalprob(gdist)
        coal_rate = np.sum(probs)
        # apply minimum nonzero coal_rate threshold
        if rmin and (coal_rate < self.err):
            coal_rate = 0.0

        # return neutral fixation rate = sum of coalescence probs
        return coal_rate

    def grid_fixrates(self, lim=False, rmin=False, fout=None):
        """get fixation rates at grid points. remove runs of 0s"""
        # gpos grid spaced by relative error parameter (start from segment end)
        gpos = self.grid_gpos[1:]
        # get approx physical grid positions and unique sort idx
        ppos, idx = np.unique(self.gmap.interp_pos(gpos), return_index=True)
        # remove gpos corresponding to ppos duplicates with unique sort idx
        gpos = gpos[idx]
        # extend physical positions to cover chromosome boundaries
        ppos[-1] = self.gmap.chlen
        # assert positive length
        assert np.all(ppos[1:] - ppos[:-1] > 0)

        # initialize containers, counters and flags
        cvals, blks = [], []  # lists for cs values and block lengths
        prev = 0  # position of previous block end
        cs_zero = False  # flag indicating previous cs score of 0
        rfix_prev = 0.0  # fixation rate (cs)

        # percentage completion printout to stderr
        pct = int(0.01 * gpos.size)  # 1% of total grid points
        stderr.write('creating sweep map: ')
        stdout.flush()

        # loop through grid sites
        for (i, (gp, pp)) in enumerate(izip(gpos, ppos)):
            # track percent progress
            if i%pct == 0:
                msg = '.' if i%(pct*25) else '{}%'.format(i/pct)
                stderr.write(msg)
                stdout.flush()

            # # if fix rate > 0 or an initial 0, create a new block
            # rfix = self.fixrate_point(gp, lim=lim, rmin=rmin)
            # if (rfix > 0) or (not cs_zero):
            #     # record cs value and block length
            #     cvals.append(rfix)
            #     blks.append(pp - prev)
            #     # flag cs value as zero/nonzero to skip consecutive 0s
            #     # cs_zero = bool(rfix)
            #     cs_zero = (rfix == 0)
            #     # update previous position
            #     prev = pp

            # # skip over consecutive 0s
            # else:
            #     continue

            # get fix rate
            rfix = self.fixrate_point(gp, lim=lim, rmin=rmin)
            # if the rate is the same as the previous rate, continue
            if rfix == rfix_prev:
                continue

            # if rates do not match, check relative difference
            else:
                # if last value was zero, create block by default
                if rfix_prev == 0:
                    blks.append(pp - prev)
                    cvals.append(rfix)
                    prev = pp
                    rfix_prev = rfix
                # otherwise, check relative change from last value
                else:
                    # if relative change > err rate, create block
                    diff = abs(rfix-rfix_prev)
                    if diff / rfix_prev > self.err:
                        blks.append(pp-prev)
                        cvals.append(rfix)
                        prev = pp
                        rfix_prev = rfix

        # add out final segment if consecutive 0s leave chrom unfinished
        end_blk = self.gmap.chlen - np.sum(blks)
        if end_blk > 0:
            cvals.append(rfix_prev)
            blks.append(end_blk)

        # completed counter message and newline
        stderr.write('\n')
        stdout.flush()

        # convert cs values and blocks to arrays, assert positve block lengths
        cvals = np.array(cvals)
        blks = np.array(blks)

        # if outfile name is given, write results to file
        if fout:
            np.savez_compressed(fout, cvals=np.column_stack((cvals, blks)))

        return cvals, blks

    def coalprob(self, r_xy):
        """coalescence prob with rec prob r_xy to sweep site """
        return np.exp(-r_xy * self.tau)

    @property
    def chrom(self):
        return self.anno.chrom

    @property
    def gmap(self):
        return self.anno.gmap

    @property
    def tau(self):
        """diffusion fixation time for a beneficial substitution"""
        ns = 4.0 * self.ne * self.coef  # 4Nes term
        ns_inv = 1.0 / ns  # 4Nes^-1
        gamma = -psi(1)  # Euler-Mascheroni constant "gamma" from psi function
        return 2.0 * (np.log(ns) + gamma - ns_inv) / self.coef

    @property
    def grid_step(self):
        """max gdist between neutral sites where fixrate changes by < 2*err"""
        return 2*self.err / self.tau

    @property
    def grid_gpos(self):
        """grid of neutral gmap positions for estimating CS within +/-err"""
        g_start, g_end = self.gmap.gmap_limits
        return np.arange(g_start, g_end, self.grid_step)


def main():
    if root_dir.startswith('/Users/davidmurphy'):
        ch = 'chr22'
        expo = -2

    else:
        if len(argv) != 3:
            print 'usage: cscalculator <chrom> <coef>'
            exit()
        ch = argv[1]
        expo = float(argv[2])

    coef = 10**expo
    an = 'nonsyn'
    tk = 'pr95.cleanrun'
    ne = 2e4
    cdir = '{}_fixed_cmap'.format(an)

    # file paths and params
    cst = ChromStruct(ch, tkn=tk, fcl=an, cdir=cdir)
    cfile = cst.cmap_file(an, coef)
    afile = cst.cs_target(anno=an).replace('.npz', '.npy')

    # annotation and map classes
    gmp = GeneticMap(ch, cst.gmap_files, gscale=0.01)
    ano = AnnoPoints(gmp, afile)

    # initialize CS calculator and save grid of fixation rates
    csc = CSCalculator(ano, coef, ne)

    start = dt.now()
    c = csc.grid_fixrates(lim=True, fout=cfile)[0]
    t = dt.now() - start
    msg = 'cmap length: {} segments\nruntime: {}\n'.format(len(c), t)
    stderr.write(msg)
    stdout.flush()

    # # FIX SAVED FILES WITH LEADING 1BP SEGMENTS:
    # cmap = np.load(cfile)['cvals']
    # assert cmap[0, 1] == 1
    # cmap = cmap[1:]  # cut off first segment
    # cmap[0, 1] += 1  # add the bp to second segment
    # np.savez_compressed(cfile, cvals=cmap)  # re-save

    # # TESTING DIFFERENT LIMIT CASES:
    # nolim = csc.grid_fixrates()
    # lim = csc.grid_fixrates(lim=True)
    # tdata = root_dir + '/lsm_python/temp_data/'
    # csc.grid_fixrates(rmin=True, fout=tdata+'rmin.npz')


if __name__ == '__main__':
    main()
