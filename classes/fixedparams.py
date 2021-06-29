from dictlike import DictLike
import numpy as np

__author__ = 'davidmurphy'


class FixedParams(DictLike):
    """
    A class that defines all of the fixed parameters
    """
    def __init__(self, bn=0, cn=0, **kwargs):
        """
        Dictlike struct of the basic fixed parameters, indices and boundaries.
        These will NOT change during the inference run or any other process.
        :param bn: number of total parameters values for bs
        :param cn: number of total parameters values for cs
        """
        super(FixedParams, self).__init__()
        # indices for BS, CS and TAU parameters in the parameters vector
        self.bn, self.cn = bn, cn
        self.bi, self.bj = (0, self.bn)
        self.ci, self.cj = (self.bn, self.bn + self.cn)
        # self.xi = self.bn + self.cn
        # self.ti = self.bn + self.cn + 1
        # self.num_params = self.bn + self.cn + 2
        self.ti = self.bn + self.cn
        self.num_params = self.bn + self.cn + 1

        # boundaries used in certain configurations of CLH
        self.min_bsx = None
        self.min_red = None  # a boundary value for the reduction
        self.min_pi0 = None  # boundary value for neutral pi
        self.max_pii = None  # a boundary value for predicted pi (pii)
        self.min_bs = None
        self.cth = None  # threshold for fraction of pii to weight with pi_avg
        # self.min_bwt = None  # weight below which to round to 0
        self.u_fix = 7.4e-8  # udel used in precalc (from McVicker result)
        self.u_max = 5e-7  # a limit for u_del (only used in some applications)
        # self.u_min = 1e-12  # a limit for udel min
        self.min_log_t = -9.999
        self.bsmax, self.bsmin = -6.0, -10.0  # max, min log(u_del) for BS
        self.csmax, self.csmin = 0.0, -10.0  # max, min log(beneficial) for CS
        self.taumax, self.taumin = -1.0, -1.0  # max min tau default
        self.bs_init = -9.0  # init BS u_del
        self.cs_init = -8.0  # init CS exponent
        self.th_init = 0.01  # initial error rate c
        self.tau_init = -1.0
        # upper/lower param bound tuples (used in some optimizations)
        bs_bounds = [(self.bsmin, self.bsmax) for _ in xrange(self.bn)]
        cs_bounds = [(self.csmin, self.csmax) for _ in xrange(self.cn)]
        tau_bounds = [(self.taumin, self.taumax)]
        self.min_max_bounds = bs_bounds + cs_bounds + tau_bounds
        # TODO: temporary!
        self.meandiv = 0.0
        # initialize from keyword args
        self.setitems(items=kwargs.items())

    @property
    def initial_parameters(self):
        """
        Follow a specific procedure to produce an array of each free param in
        the run:
        1) check the first BS annotation in self.bs_annos
        2) add a free param for each selection coefficient given for this anno
        3) go to the next BS annotation and repeat (2) until BS annos are
        complete
        4) repeat steps (1-3) for CS annotations
        5) add a single param for tau
        Example:
        - An inference has 3 BS annos and 2 CS annos
        - All BS annos have 6 selection coefficients and and all CS annos have
        4 selection coefficients
        - The inference has 3x6 + 2x4 + 1 = 27 free parameters
        NOTE: Initial param values are the 'init' values for bs, cs, and tau
        designated in FixedParams
        :return: the initial vector of params using init values
        """
        parameters = np.zeros(shape=self.num_params, dtype='f8')
        parameters[self.bi:self.bj] = self.bs_init
        parameters[self.ci:self.cj] = self.cs_init
        # parameters[self.xi] = self.th_init
        parameters[self.ti] = self.tau_init  # note: update at pre-process step
        return parameters
