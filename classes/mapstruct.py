from itertools import izip
from datetime import datetime
from sys import stdout, stderr
from multiprocessing import Pool
from datasampler import DataSampler
from functions import weighted_variance
from functions import time2str, str2time
from scipy.optimize import minimize as minimize
from runstruct import RunStruct, ChromStruct, np, os


# from memory_profiler import profile

__author__ = 'davidmurphy'


def _trynpy(fname):
    """try loading numpy file, print detailed message if memory overflow"""
    try:
        return np.load(fname)
    except MemoryError:
        message = 'memory overflow while loading {}\n'.format(fname)
        stderr.write(message)
        stdout.flush()
        exit(1)


def _kwdpcalc(pdict):
    """wrapper function that calls self.calc on a key=arg dict"""
    return _pcalc(**pdict)


def _pcalc(prm=None, fx=None, hom=None, het=None, u=None, bs=None, cs=None,
           s=None, **kwargs):
    """Calculate log CLH with direct inputs"""

    log_lh = 0.0
    # message = 'pcalc on data index={}'.format(i)
    # stderr.write(message)
    # stdout.flush()

    # boundary conditions:
    min_red = fx.min_red
    min_pi0 = fx.min_pi0
    max_pii = fx.max_pii

    for _ in xrange(1):

        # reduction in pi from bs -> exp(dot(bs, bs_weights))
        if bs is not None:
            # convert bs params to udel/site for each bmap
            uvec = np.power(10, prm[fx.bi:fx.bj])
            bwt = uvec / fx.u_fix
            bsx = np.exp(np.dot(bs, bwt))
        else:
            bsx = np.ones(shape=len(u))

        # reduction in pi from cs -> dot(cs, cs_params)
        if cs is not None:
            cwt = np.power(10, prm[fx.ci:fx.cj])
            csx = np.dot(cs, cwt)
        else:
            csx = np.zeros(shape=len(u))

        # theta is initially mean pi scaled by mutation rate variation estimate
        theta = u / prm[fx.ti]
        # use meandiv instead of mean(u) -> more precise
        mean_theta = fx.meandiv / prm[fx.ti]
        # free param pi0 = expected het w/o selection
        pi0 = np.maximum(theta / (1.0 + theta), min_pi0)
        # red = overall reduction in pi from bs & cs:
        red = (1 + mean_theta) / (1 / bsx + csx + mean_theta)
        # pii = pi after reduction from BS & CS
        pii = np.minimum(pi0 * np.maximum(red, min_red), max_pii)

        # standard logLH calculation for the segment of neutral sites
        if s is None:
            log_lh = hom * np.log(1 - pii) + het * np.log(pii)
        else:
            log_lh = s * hom * np.log(1 - pii) + s * het * np.log(pii)

    # calculate the scaled -sum in advance so only one term sent back to main
    sum_llh = -1e5 * np.sum(log_lh)

    return sum_llh


class CompositeLikelihood(RunStruct):
    pass


class MapStruct(ChromStruct):
    """
    A data structure inheriting the properties of a ChromStruct, although in
    practice the MapStruct combines data from across the autosomes.
    """
    # @profile
    def __init__(self, seg=False, div=False, pool=False, **kwargs):
        """
        Initialize ChromStruct parent with kwargs and fixed chromosome
        :param seg: flag that indicates that segments should be loaded
        :param div: flag that indicates that divergence data should be loaded
        :param multp: flag that indicates multi-threading to be used
        """
        # init parent class
        super(MapStruct, self).__init__(chrom='chr1', **kwargs)

        # # DEBUG
        # message = 'initializing MapStruct\n'
        # stderr.write(message)
        # stdout.flush()

        # neutral data for each chrom as a list and save concatenated arrays
        ntarrs = self._getmaps(self.nt_files, mask=False, aslist=True)

        # # DEBUG
        # message = 'neutral arrays loaded into list\n'
        # stderr.write(message)
        # stdout.flush()

        # create "has-data" mask for each neutral data array
        self.msk = [(a[:, 0] > 0) for a in ntarrs]

        # # DEBUG
        # message = 'neutral array masks created\n'
        # stderr.write(message)
        # stdout.flush()

        # save unmasked array lengths
        self.lens = [len(a) for a in ntarrs]
        # use lengths to calculate boundaries of unmasked concatenated array
        self.bnds = [0] + [s for s in np.cumsum(self.lens)]
        # record the number of segments in the data
        self.cnt = sum(self.lens)

        # get lengths, boundaries and total seg count for masked arrays as well
        self.mlens = [np.sum(m) for m in self.msk]
        self.mbnds = [0] + [s for s in np.cumsum(self.mlens)]
        self.mcnt = sum(self.mlens)

        # concatenate mask across chroms
        self.ms = np.concatenate(self.msk)
        # reduce the final neutral data array to segments containing data only
        self.hom, self.het = np.concatenate(ntarrs)[self.ms].T

        # # DEBUG
        # message = 'concatenated mask, concatenated neutral array created\n'
        # stderr.write(message)
        # stdout.flush()

        # load bs & cs arrays if they exist
        self.bs = self._logb if self.bnum else None
        self.cs = self._getmaps(self.cs_files) if self.cnum else None

        # # DEBUG
        # message = 'BS/CS maps arrays loaded\n'
        # stderr.write(message)
        # stdout.flush()

        # load umap
        self.u = self._getmaps(self.u_files)

        # # DEBUG
        # message = 'concatenated u vector created\n'
        # stderr.write(message)
        # stdout.flush()

        # load segments and divergence data optionally
        if seg:
            self.seg = self._getmaps(self.seg_files, mask=False, aslist=True)
        else:
            self.seg = None
        if div:
            self.div = self._getmaps(self.ndiv_files)
        else:
            self.div = None

        # bootstrap/jackknife mask if applicable
        # TODO: turn DataSampler into a ChromStruct
        if self.vars.sample_type is not None:
            self.s = DataSampler(clone=self).get_mask()[self.ms]
        else:
            self.s = None

        # initialize pool of processors for multi-threading
        if pool:
            # self.pool = 1
            self.pool = Pool(processes=self.vars.num_cores)
            self._inputs = self.pcalc_inputs(n_pieces=self.vars.num_cores)

            # # DEBUG
            # message = 'pcalc_inputs created\n'
            # stderr.write(message)
            # stdout.flush()

        else:
            self.pool = None
            self._inputs = None

        # # list of optimization log output that dumps every 200 lines
        # self.cached_params = []

    """
    LIKELIHOOD CALCULATION AND PREDICTED DIVERSITY MAP
    --------------------------------------------------
    """
    def optimize(self):
        """
        Notes on some of the optimization methods used:
        -- Nelder-Mead is fast and takes large steps but it cannot be bounded
        so it may return negative weights
        -- In practice we fix the negative weights in the log_likelihood
        function and call them either 0 or min_w
        -- Truncated Newton is close to the 'active-set' option used in MATLAB
        """
        # # delete init file if it exists, save new copy with timestamp
        # if self.init_file is not None:
        #     os.remove(self.init_file)

        # record run start time format=YYMMDDHHMMSS
        self.record_timestamp()
        # affix start time to init_file to
        self.save()

        # switch off the 'complete' flag before entering optimization
        self.vars.complete = False

        # loop through optimization methods
        for optimizer_id in self.optimizers:

            # # DEBUG
            # message = 'begin {} minimization\n'.format(optimizer_id)
            # stderr.write(message)
            # stdout.flush()

            # store start/end times for each optimization method
            start_time = datetime.now()

            # get optimizer data structure for 'optimizer_id' to parameterize
            optimizer = self[optimizer_id]

            # run optimizers in sequence and use each optimizer result as
            # the set of initial params for the following optimizer
            optimization = minimize(fun=self.log_likelihood, x0=self.params,
                                    **optimizer.methodargs)

            # process OptimizeResult object
            optimizer.optimizer_results(result=optimization)

            # get final params, use as x0 for next method
            self.params = optimizer.x

            # update best method, likelihood and param values
            if optimizer.fun < self.stat.best_lh:
                self.stat.best_lh = optimizer.fun
                self.stat.best_method = optimizer.method
                self.stat.best_params = optimizer.x

            # save the method run time and cumulative run time
            optimizer.runtime = time2str(datetime.now() - start_time)
            hours = str2time(optimizer.runtime).total_seconds() / 3600.0
            self.stat.total_time += round(hours, 2)

        # write the final output to the optimization log
        if len(self.cached_params) > 0:
            self.log_params()

        # mark the run as complete
        self.vars.complete = True

        # calc final stats on predicted maps and params, save to final_file
        self.final_stats()
        self.save(txt_file=self.final_file)

        # delete timestamped initialization file
        # os.remove(self.txt_file)

    def log_likelihood(self, params):
        """
        Calculates the likelihood of the neutral data given a weighted
        combination of the B-maps
        :param params: a variable copy of free params attribute that is
        manipulated by the optimizer and used to re-set the internal attribute
        param values
        :return neg_log_p, dneg_log_p: the negative log probability and vector
        of 1st derivatives for free params
        """
        # if self.stat.function_calls > 20:
        #     exit(1)

        # increment the function call counter
        self.stat.function_calls += 1

        # # DEBUG
        # tests = len(self.params), len(params)
        # message = 'len(self.params)={} len(params)={}\n'.format(*tests)
        # stderr.write(message)
        # stdout.flush()

        # update internal params
        self.params = params

        # DEBUG
        st = datetime.now()
        message = 'starting {} LLH calc... '.format(self.process_type)
        stderr.write(message)
        stdout.flush()

        # regular single thread function
        if self.pool is None:
            # logLH and first derivatives for free parameters
            lh_vector, dlh_vector = self.calc(**self.dict)
            # multiply sum by 1e5 to avoid underflow
            nlp = -1e5 * np.sum(lh_vector) / self.stat.pairs

        else:
            # update params in each subset of data
            self._updateprm()
            # sum the piecemeal (-1e5)-scaled logLH sums from multi-thread
            pcalc_outputs = self.pool.map(_kwdpcalc, self._inputs)
            nlp = sum(pcalc_outputs) / self.stat.pairs

        # DEBUG
        en = datetime.now()
        message = 'LLH={}, time={}\n'.format(nlp, en-st)
        stderr.write(message)
        stdout.flush()

        # record params, -sum(logCLH) and udel for each bs anno
        self.cur_params(nlp)

        return nlp

    def calc(self, hom=None, het=None, u=None, bs=None, cs=None, s=None,
             **kwargs):
        """Calculate log CLH and first derivatives"""
        if (bs is None) and (cs is None):
            raise TypeError('calc() requires bs or cs (os both)')

        # boundary conditions:
        min_red = self.fixed.min_red
        min_pi0 = self.fixed.min_pi0
        max_pii = self.fixed.max_pii

        # calculate effects of BS & CS
        bs = self.bs_prediction
        cs = self.cs_prediction

        # theta is initially mean pi scaled by mutation rate variation estimate
        theta = u / self.params[self.fixed.ti]

        # use meandiv instead of mean(u) -> more precise
        mean_theta = self.stat.meandiv / self.params[self.fixed.ti]
        # free param pi0 = expected het w/o selection

        # NOTE: "pi0", "pii" and "red" are given boundary conditions to
        # avoid inference crash for some optimization algorithms used in the
        # MATLAB version (Elyashiv et al., 2016) --> are these still needed and
        # is it possible that they could mess up the current optimization?
        pi0 = np.maximum(theta / (1.0 + theta), min_pi0)

        """
        red = overall reduction in pi from bs & cs:
        
                       (1 + theta)
                   ___________________    Eq. 1 (Elyashiv et al., 2016)
                
                   (1/bs + cs + theta)
        
        """
        red = (1 + mean_theta) / (1/bs + cs + mean_theta)

        # pii = diversity after reduction from background selection and
        # selective sweeps.
        pii = np.minimum(pi0 * np.maximum(red, min_red), max_pii)

        # just return the predicted pi vector if results mode is set in mst
        if self.vars.complete:
            return pii, pi0

        # standard logLH calculation for the segment of neutral siteself.s
        if s is None:
            log_lh = hom * np.log(1 - pii) + het * np.log(pii)
        else:
            log_lh = s * hom * np.log(1 - pii) + s * het * np.log(pii)

        return log_lh, []

    """
    STATS ON DATA AND FINAL MAP
    ---------------------------
    """
    def initial_stats(self):
        """calculate basic statistics from the data"""

        # load u, segments, divergence maps
        self.u = self._getmaps(self.u_files)
        self.seg = self._getmaps(self.seg_files, mask=False, aslist=True)
        self.div = self._getmaps(self.ndiv_files)

        # record descriptive stats about segments used in the infererence
        self.stat.total_segs = self.cnt
        self.stat.used_segs = self.mcnt
        bases_covered = sum(s[m].sum() for (s, m) in zip(self.seg, self.msk))
        self.stat.covered_bases = bases_covered

        # create a vector from the bootstrap if applicable otherwise use ones
        if self.s is not None:
            bts = self.s
        else:
            bts = np.ones(shape=self.mcnt, dtype='f8')

        # count up mono/poly site totals
        poly = []
        for ch in self.chroms:
            fpoly = self.neut_poly.replace(self.chrom, ch)
            poly.append(np.load(fpoly)['neutpoly'])
        self.stat.indv = max(np.max(np.sum(ar[:, 1:], axis=1)) for ar in poly)
        self.stat.sites = self.nsites.sum()
        self.stat.poly = sum(len(ar) for ar in poly)
        self.stat.mono = self.stat.sites - self.stat.poly

        # diversity from polymorphism
        self.stat.homs = np.sum(self.hom * bts)
        self.stat.hets = np.sum(self.het * bts)
        self.stat.pairs = self.stat.homs + self.stat.hets
        self.stat.meanpi = self.stat.hets / self.stat.pairs  # genome-wide pi
        unweighted = self.het / np.sum(self.hom+self.het)
        self.stat.varpi = weighted_variance(x=unweighted, w=bts)

        # weight mean and var by sites/segment and bootstrap for each segment
        wts = bts * self.nsites

        # mean divergence by substitution count per neutral site
        divx = self.div[:, 1] / self.nsites
        self.stat.meandiv = np.average(divx, weights=wts)
        self.stat.vardiv = weighted_variance(x=divx, w=bts)
        # mean u from u estimates
        self.stat.meanu = np.average(self.u, weights=wts)
        self.stat.varu = weighted_variance(x=self.u, w=wts)

        # count up the number of sites covered by bs targets
        self.stat.bsanno_cov = []  # init empty for safety
        for a in self.bs_annos:
            asum = 0
            for ch in self.chroms:
                fname = self.bs_targets[a].replace(self.chrom, ch)
                start, end = np.loadtxt(fname, usecols=(1, 2), dtype='u8').T
                asum += np.sum(end - start)
            self.stat.bsanno_cov.append((a, asum))

        # set initial tau and boundaries
        self.fixed.tau_init = self.stat.meandiv / self.stat.meanpi
        self.fixed.taumin = 0.1 * self.fixed.tau_init
        self.fixed.taumax = 10.0 * self.fixed.tau_init

        # set tau param ONLY if it is set to the default initial value (-1.0)
        if self.params[self.fixed.ti] == -1.0:
            self.params[self.fixed.ti] = self.fixed.tau_init

        # record the initial logLH before optimization
        self.vars.complete = False
        self.stat.initial_lh = self.log_likelihood(params=self.params)

        # re-write init file with new stats included
        self.save()

    def final_stats(self):
        """calculate stats from inferred params and prediction maps"""
        # run must be completed
        assert self.vars.complete

        # create a vector from the bootstrap if applicable otherwise use ones
        if self.s is not None:
            bts = self.s
        else:
            bts = np.ones(shape=self.mcnt, dtype='f8')

        # weight segments by number of sites times bootstrap coefficient
        wts = bts * self.nsites
        # get predicted map of diversity and expected map under neutrality
        assert self.vars.complete
        pred, pi0 = self.calc(**self.dict)

        # use weights for summary stats
        self.stat.meanpi0 = np.average(pi0, weights=wts)
        self.stat.varpi0 = weighted_variance(x=pi0, w=wts)
        self.stat.meanpred = np.average(pred, weights=wts)
        self.stat.varpred = weighted_variance(x=pred, w=wts)
        # calculate overall average reduction in diversity
        self.stat.diversity_reduction = 1 - self.stat.meanpi/self.stat.meanpi0

        # clear existing results so there is no overwrite
        self.stat.uvec = []
        self.stat.utot = []
        self.stat.upmf = []
        self.stat.urel = []
        self.stat.tmean = []
        self.stat.tstd = []

        # BS param summaries
        i = 0
        for t_vec in self.bdfe:
            j = i + len(t_vec)
            u_vec = 10 ** self.params[i:j]
            u_tot = np.sum(u_vec)
            u_pmf = u_vec / u_tot
            u_rel = u_vec / self.fixed['u_fix']
            t_mean = np.average(t_vec, weights=u_pmf)
            t_var = np.average(t_vec ** 2, weights=u_pmf) - t_mean ** 2
            self.stat.uvec.append(u_vec)
            self.stat.utot.append(u_tot)
            self.stat.upmf.append(u_pmf)
            self.stat.urel.append(u_rel)
            self.stat.tmean.append(t_mean)
            self.stat.tstd.append(np.sqrt(t_var))
            i = j

        self.save()

    def prediction(self, umap=None):
        """return model predictions for umap and current params"""
        # set "complete" flag to True so calc returns pii and pi0
        self.vars.complete = True
        kwargs = self.dict
        if umap is not None:
            kwargs['u'] = umap
        return self.calc(**kwargs)[0]

    def window_index(self, window):
        """return indices dividing data into a grid of window-size blocks"""
        if self.seg is None:
            self.seg = self._getmaps(self.seg_files, mask=False, aslist=True)
        # get the cumulative positions across autosomes
        cumpos = np.cumsum(np.concatenate(self.seg))[self.ms]
        # the grid within which to partition subsets of the data
        grid = np.arange(window, cumpos[-1] + window, window)
        # return indices that sort data into the grid
        return np.searchsorted(a=cumpos, v=grid)

    def split_index(self, n_pieces, masked=True):
        """
        A function that returns a set of start:end indices that will split
        concatenated data arrays into "n_pieces" different parts, e.g., if
        self.bs has len=1000, split_index(n_pieces=4) would return the indices
        [(0,250), (250, 500), (500, 750), (750, 1000)]
        :param n_pieces: number of subsets to break the data into
        :param masked: if "masked" flag is raised, use the number of segments
        AFTER masking (self.mcnt) as the basis for the division; otherwise
        use the total number of segments (self.cnt) before masking as the
        basis for the division
        :return:
        """
        # set the number of segments used as basis for sub-division
        nsegs = self.mcnt if masked else self.cnt

        # get the step-size to evenly subdivide concatenated data
        step = int(nsegs / n_pieces)

        # the there should be fewer than "n_pieces" segments cut off
        assert int(nsegs % n_pieces) < n_pieces

        # use cumulative sum of steps to get the end indices
        jdx = np.cumsum(np.full(n_pieces, step, dtype=int))
        idx = jdx - step  # start indices are just (upper indices - step)
        jdx[-1] = nsegs  # adjust for segments cutoff by subdivision

        # join start and end indices into (n_pieces, 2) array
        sidx = np.column_stack((idx, jdx))

        # return (n_pieces x 2) array of [start, end] indices
        return sidx

    def _getmaps(self, fname, mask=True, aslist=False):
        """list or concatenate arrays for fname for all chroms"""
        arrs = []
        # mask out rows where there is no neutral diversity data
        if mask:
            assert len(self.msk) == len(self.chroms)
            for (ms, ch) in izip(self.msk, self.chroms):
                f = fname.replace(self.chrom, ch)
                ar = _trynpy(f)[ms].astype('f8')
                arrs.append(ar)

                # # DEBUG
                # message = '{} loaded and masked\n'.format(f)
                # stderr.write(message)
                # stdout.flush()

        # keep all rows
        else:
            for ch in self.chroms:
                f = fname.replace(self.chrom, ch)
                ar = _trynpy(f).astype('f8')
                arrs.append(ar)

                # # DEBUG
                # message = '{} loaded\n'.format(f)
                # stderr.write(message)
                # stdout.flush()

        # return a list with one array per chromosome
        if aslist:
            return arrs
        # concatenate arrays across chromosomes
        else:
            return np.concatenate(arrs)

    @property
    def bs_prediction(self):
        """return the predicted effects of BS for the current parameters"""
        if self.bnum:
            # convert bs params to udel/site for each bmap
            urates = np.power(10, self.params[self.fixed.bi:self.fixed.bj])
            weights = urates / self.fixed.u_fix
            bs = np.exp(np.dot(self.bs, weights))
            return bs
        else:
            return np.ones(shape=self.mcnt)

    @property
    def cs_prediction(self):
        """return the predicted effects of CS for the current parameters"""
        if self.cnum:
            cwt = np.power(10, self.params[self.fixed.ci:self.fixed.cj])
            cs = np.dot(self.cs, cwt)
            return cs
        else:
            return np.zeros(shape=self.mcnt)

    @property
    def _logb(self):
        """take log of bmaps rescaled to 0.0-1.0"""
        bmap_array = np.maximum(1, self._getmaps(self.bs_files))
        rescaled_bmap = bmap_array / float(self.bscl)
        return np.log(rescaled_bmap)

    def pcalc_inputs(self, n_pieces):
        """stratify data to distribution CLH calculation across cores"""
        # # DEBUG
        # message = 'stratify data into {} bins\n'.format(self.vars.num_cores)
        # stderr.write(message)
        # stdout.flush()

        # get sub-division indices
        idx = self.split_index(n_pieces)

        # create sub-arrays from data matrices using splitting indices
        pcalc_inputs = []
        for (i, j) in idx:
            # required poly & div matrices
            hom, het, u = self.hom[i:j], self.het[i:j], self.u[i:j]

            # optional bs, cs, s matrices
            bs = self.bs[i:j] if self.bs is not None else None
            cs = self.cs[i:j] if self.cs is not None else None
            s = self.s[i:j] if self.s is not None else None

            # store constant/stratified inputs together as a tuple
            pcalc_inputs.append([hom, het, u, bs, cs, s])
            # pdict = dict(hom=hom, het=het, u=u,
            #              bs=bs, cs=cs, s=s)
            # pcalc_inputs.append(pdict)

        return pcalc_inputs

    def _updateprm(self):
        """update the params in each split data dict"""
        for pdict in self.pcalc_inputs:
            pdict['prm'] = self.params

    @property
    def chrom_tag(self):
        """return vector of chrom index for each row of nt"""
        enum = enumerate(self.mlens, start=1)
        return np.concatenate([np.full(l, i, dtype='u1') for (i, l) in enum])

    @property
    def process_type(self):
        return 'single-thread' if self.pool is None else 'multi-thread'

    """
    POSITIONS, RANGES, SITE COUNTS
    ------------------------------
    """
    @property
    def nt(self):
        """stack hom, het columns and return array"""
        return np.column_stack((self.hom, self.het))

    @property
    def pos(self):
        """seg lengths converted to cum pos on each autosome separately"""
        if self.seg is None:
            self.seg = self._getmaps(self.seg_files, mask=False, aslist=True)
        pos_list = [np.cumsum(s, dtype='u4') for s in self.seg]
        return np.concatenate(pos_list)[self.ms]

    @property
    def ranges(self):
        """the ranges of each data segment in chromosome bed coordinates"""
        if self.seg is None:
            self.seg = self._getmaps(self.seg_files, mask=False, aslist=True)
        start_list = self.pos - np.concatenate(self.seg)[self.ms]
        # stack lower and upper positions from each segment (masked)
        return np.column_stack((start_list, self.pos))

    @property
    def nsites(self):
        """number of neutral sites per segment across neutral segments"""
        if self.div is None:
            self.div = self._getmaps(self.ndiv_files)
        return self.div[:, 0]

    @property
    def cumsites(self):
        """cumulative neutral sites/segment for each chrom"""
        # lower, upper indices partitioning masked segments into chroms
        idx = izip(self.mbnds[:-1], self.mbnds[1:])
        count_list = [np.cumsum(self.nsites[i:j]) for (i, j) in idx]
        return np.concatenate(count_list)

    @property
    def slices(self):
        """indices partitioning individual sites into segments across chroms"""
        return np.column_stack([self.cumsites - self.nsites, self.cumsites])

    """
    GENETIC MAP QUANTITIES
    ----------------------
    """
    @property
    def cmsegs(self):
        """return the distance in cM between edges of each segment"""
        # all ranges masked
        ranges = self.ranges
        cms = []
        for (i, ch) in enumerate(self.chroms):
            # load map for each chrom
            fname = self.gmap_files.replace(self.chrom, ch)
            pos, gpos = np.loadtxt(fname, usecols=(0, 2)).T
            # interpolate cm values for segment edges
            edges = ranges[self.mbnds[i]:self.mbnds[i+1]]
            cmi, cmj = np.interp(edges, pos, gpos).T
            # record the segment length in cM
            cms.append(cmj - cmi)

        return np.concatenate(cms)

    @property
    def cmmbsegs(self):
        """return the cM/Mb between each segment"""
        return 1e6 * self.cmsegs / self._getmaps(self.seg_files)

    @property
    def uconst(self):
        """make u a vector constant equal to meandiv"""
        return np.full(self.mcnt, self.stat.meandiv)

    """
    PREDICTION SORTING INDICES
    --------------------------
    """
    @property
    def prediction_sorting_index(self):
        """sorting index orders segments by diversity prediction (low-high)"""
        return self.prediction(umap=self.uconst).argsort()

    # def get_bmap(self, chrom, cum=False):
    #     """
    #     Get the completed bmap (mcvicker style) from a completed experiment
    #     :param chrom: chromosome to map
    #     :param cum: flag to return the cumulative positions rather than
        # segment lengths for bs values
    #     :return bmap: [seg, b] formatted array
    #     """
    #     assert isinstance(self, RunStruct)
    #     idx = self.chroms.index(chrom)
    #     barr = np.log(np.maximum(np.load(self.bs_files[idx]), 1) / 100.0)
    #     sarr = np.load(self.seg_files[idx])
    #     if cum:
    #         sarr = np.cumsum(sarr)
    #     uvec = np.power(10, self.params[self.fixed.bi:self.fixed.bj])
    #     bwt = uvec / self.fixed.u_fix
    #     bs = np.exp(np.dot(barr, bwt))
    #     bs *= 100.0
    #
    #     return np.column_stack([bs, sarr])

