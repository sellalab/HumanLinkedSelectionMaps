from dictlike import DictLike

__author__ = 'davidmurphy'


class Variables(DictLike):
    """
    A DictLike to manage variables from the inference like bootstrap_index,
    optimizer methd:name pairs, etc.
    """
    def __init__(self, **kwargs):
        """Load empty Variables class and set kwargs attributes (if any)"""
        super(Variables, self).__init__()

        # number of cores to use for multiprocessing
        self.num_cores = 1

        # a flag to calculate the derivative for the CLH function
        self.calc_deriv = False

        # phastCons min, max score range (default is [0, 0]) for NEUTRAL sites
        self.pmin = 0.0
        self.pmax = 0.0
        self.percentile = 95

        # "PASS" symbol for 1000 Genomes sequencing call mask
        self.callmask_key = 'P'

        # optional filters
        self.filter_cpg = False
        self.filter_con = True

        # params for estimating mutation rate variation
        self.muest_window = 2e4  # bp-window for rate mutation estimates
        self.min_bases = 100  # min bases/window for mutation rate estimate

        # flag to use a constant mutation rate instead of divergence estimates
        self.mu_const = False

        # params for collated plots
        self.collated_bin = 1.25e-5  # cM width of bins in collated plots
        self.collated_span = 0.4  # max(abs(focal-neutral)) for collated plots

        # sampling strategy (e.g., bootstrap or jackknife)
        self.sample_type = None

        # down sampling flag
        self.down_sample = False
        self.down_sample_pct = 0.01

        # bootstrap params
        self.bootstrap_window = 1e6
        self.bootstrap_samples = 0
        self.bootstrap_index = 0

        # jackknife params
        self.use_jackknife = False
        self.jackknife_window = 2e6
        self.jackknife_samples = 0
        self.jackknife_index = 0

        # simulation input params
        self.simulation_params = []

        # flags reflecting the status of the RunStruct and whether
        self.complete = False  # set to true after optimization is completed

        # initialize from keyword args
        self.setitems(items=kwargs.items())

    @property
    def multi_thread(self):
        return self.num_cores > 1
