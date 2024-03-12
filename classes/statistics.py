from data_processing.data_tools import str2time
from classes.dictlike import DictLike
import numpy as np

__author__ = 'davidmurphy'


class Statistics(DictLike):
    """
    A data structure of descriptive and analytic stats and inference results
    with basic calculator functions.
    """
    def __init__(self, **kwargs):
        """Load empty Statistics class and set kwargs attributes (if any)"""
        super(Statistics, self).__init__()
        # descriptive statistics for diversity dataset
        self.indv = 0.0  # number of individuals in the polymorphic data
        self.poly = 0.0  # polymorphic sites
        self.mono = 0.0  # monomorphic sites
        self.sites = 0.0  # total sites (mono+poly)
        self.hets = 0.0
        self.homs = 0.0
        self.pairs = 0.0
        # descriptive statistics for annotations and pre-calculated maps
        self.bsanno_cov = []
        self.csanno_cnt = []
        self.total_segs = 0.0  # all segments after compression
        self.used_segs = 0.0  # segments after masking segs w/o data
        self.bases_covered = 0.0  # total length of used segs
        # diversity and divergence statistics from data
        self.meanpi = 0.0
        self.varpi = 0.0
        self.meandiv = 0.0  # mean divergence from total subs / total sites
        self.vardiv = 0.0
        self.meanu = 0.0  # mean divergence from mutation rate estimates
        self.varu = 0.0
        # statistics on diversity predictions
        self.meanpred = 0.0
        self.varpred = 0.0
        self.meanpi0 = 0.0
        self.varpi0 = 0.0
        # average diversity reduction due to selection
        self.diversity_reduction = 0.0
        # inferred params (e.g., DFE), param summaries (e.g., total udel)
        self.uvec = []
        self.utot = []
        self.upmf = []
        self.urel = []
        self.tmean = []
        self.tstd = []
        self.avec = []
        self.atot = []
        self.apmf = []
        self.smean = []
        self.sstd = []
        # segment counts for each chromosome: total and with masking
        self.totsegs = []
        self.msksegs = []
        # initial parameters and their CLLH
        self.init_params = None
        self.init_clh= None
        # best results out of all optimizations
        self.best_lh = 1e6
        self.best_method = None
        self.best_params = None
        # params used for simulated data and CLLH under the true params
        self.true_params = None
        self.true_clh = None
        # some run stats
        self.total_time = 0.0  # time between call to 'optimize' and its return
        self.function_calls = 0  # += 1 for every call to log_likelihood
        # R^2 values for resulting predictions across different spatial scales
        self.rsq_values = []
        # initialize from keyword args
        self.setitems(items=kwargs.items())

    def calc_stats(self, rst):
        """
        Calculate final statistics after optimization concludes
        :param rst: RunStruct controlling the concluded run
        :type rst: RunStruct
        """
        # total runtime in hours (from seconds)
        self.total_time = 0.0
        for opname in rst.optimizers:
            op = rst[opname]
            self.total_time += (str2time(op.runtime).total_seconds() / 3600.0)

        # clear each list or stats will pile up from multiple calls
        self.uvec = []
        self.utot = []
        self.upmf = []
        self.urel = []
        self.tmean = []
        self.tstd = []
        self.avec = []
        self.atot = []
        self.apmf = []
        self.smean = []
        self.sstd = []

        # BS param summaries
        i = 0
        for t_vec in rst.bdfe:
            j = i + len(t_vec)
            u_vec = 10 ** rst.params[i:j]
            u_tot = np.sum(u_vec)
            u_pmf = u_vec / u_tot
            u_rel = u_vec / rst.fixed['u_fix']
            t_mean = np.average(t_vec, weights=u_pmf)
            t_var = np.average(t_vec ** 2, weights=u_pmf) - t_mean ** 2
            self.uvec.append(u_vec)
            self.utot.append(u_tot)
            self.upmf.append(u_pmf)
            self.urel.append(u_rel)
            self.tmean.append(t_mean)
            self.tstd.append(np.sqrt(t_var))
            i = j

        # CS param summaries
        for s_vec in rst.cdfe:
            j = i + len(s_vec)
            a_vec = 10 ** rst.params[i:j]
            a_tot = np.sum(a_vec)
            a_pmf = a_vec / a_tot
            s_mean = np.average(s_vec, weights=a_pmf)
            s_var = np.average(s_vec ** 2, weights=a_pmf) - s_mean ** 2
            self.avec.append(a_vec)
            self.atot.append(a_tot)
            self.apmf.append(a_pmf)
            self.smean.append(s_mean)
            self.sstd.append(np.sqrt(s_var))
            i = j

        return None
