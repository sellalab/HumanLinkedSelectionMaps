__author__ = 'davidmurphy'

import os
from classes.runstruct import ChromStruct, root_dir, np

# copy inferred params from new map optimization results
fdir = '{}/result/final_files/sims'.format(root_dir)
fname = '00-YRI.pr95.clustinf.initial.BS1.6.CS0.0.180408214759.final.txt'
finit = '{}/{}'.format(fdir, fname)

# original run cst
old_cst = ChromStruct(chrom='chr1', init=finit)
old_cst.stat.calc_stats(old_cst)

# params from original run
mean_u = old_cst.stat.utot[0]
mean_div = old_cst.stat.meandiv
init_tau = old_cst.fixed.tau_init
true_tau = old_cst.params[-1]

"""
OPTIMIZER METHODS
"""
# optimization methods
grad_methods = ('Newton-CG', 'CG', 'BFGS')
optim_meth = ('Nelder-Mead', 'BFGS')
# optim_meth = ('BFGS',)
# optim_meth = ('Powell',)
# optim_meth = ('Newton-CG',)


class ParamGenerator(object):
    """an object that generates inference parameters"""
    def __init__(self, nba, nca, nprm, itau, ftau, udel, alpha):
        # number of bs and cs annotations
        self.nba = nba
        self.nca = nca
        # # number of bs and cs params needed for each anno
        # self.nbs = nbs
        # self.ncs = ncs
        # number of params per anno
        self.nprm = nprm
        # initial and final tau values to use for init & true params
        self.itau = itau
        self.ftau = ftau
        # total deleterious u rate (udel) and fraction beneficial subs (alpha)
        self.udel = udel
        self.alpha = alpha

    def param_sets(self, patterns=('ufm_0'), scales=(1.0,)):
        pass

    def uniform(self, pattern, scale=1.0, tau=None):
        """generates patterns of uniform parameters at different scales"""
        # set tau externally or use preset initial value
        if tau is None:
            tau = self.itau
        # don't allow alpha above 1.0 due to scaling
        cmax = np.minimum(self.alpha*scale*self.npinv, 1.0)
        # variable for half of params used in some distributions
        p_half = int(self.nprm / 2)
        # bs and cs values for half param distributions
        b_half = np.log10(self.udel * scale * self.npinv * 2)
        c_half = np.log10(cmax*2)
        # pattern 1: uniform udel/alpha across params
        if pattern == 'ufm_0':
            bprm = [np.log10(self.udel*scale*self.npinv)] * self.nprm
            cprm = [np.log10(cmax)] * self.nprm
            prms = np.array(bprm*self.nba + cprm*self.nca + [tau])
        # pattern 2: uniform udel/alpha across lower half of BS & CS params
        elif pattern == 'ufm_1':
            bprm = [b_half] * p_half + [-12.0] * p_half
            cprm = [c_half] * p_half + [-12.0] * p_half
            prms = np.array(bprm*self.nba + cprm*self.nca + [tau])
        # pattern 2: uniform udel/alpha across upper half of BS & CS params
        elif pattern == 'ufm_2':
            bprm = [-12.0] * p_half + [b_half] * p_half
            cprm = [-12.0] * p_half + [c_half] * p_half
            prms = np.array(bprm*self.nba + cprm*self.nca + [tau])
        # pattern 3: alternating uniform starting from lower half
        elif pattern == 'ufm_3':
            bprm = [b_half if (i+1)%2 else -12.0 for i in xrange(self.nprm)]
            cprm = [c_half if (i+1)%2 else -12.0 for i in xrange(self.nprm)]
            prms = np.array(bprm*self.nba + cprm*self.nca + [tau])
        elif pattern == 'ufm_4':
            bprm = [b_half if i%2 else -12.0 for i in xrange(self.nprm)]
            cprm = [c_half if i%2 else -12.0 for i in xrange(self.nprm)]
            prms = np.array(bprm*self.nba + cprm*self.nca + [tau])
        else:
            raise NameError('{} pattern does not exist'.format(pattern))

        return prms

    def random_set(self, subset, scales, resample=False, tau=None):
        """get a subset of random params for a set of scales"""
        # set tau
        if tau is None:
            tau = self.itau
        # subset indices
        i, j = subset
        rprm_set = []
        # for each scale collect the scaled subset of random params
        for scl in scales:
            rprm = self.random(scale=scl, resample=resample, tau=tau)[i:j]
            rprm_set.append(rprm)
        # combine each scaled subset into single array
        rprm_set = np.concatenate(rprm_set)

        return rprm_set

    def random(self, scale=1.0, resample=False, tau=None):
        """generates randomly weighted parameters at different scales"""
        # set tau externally or use preset initial value. create tau vector
        if tau is None:
            tau = self.itau
        tprm = np.full(25, tau)
        # divide random weights into bs and cs weights
        rwts = self.random_weights(resample=resample)
        bwts = rwts[:, :self.bi]
        cwts = rwts[:, self.bi:self.ci]
        # generate bs params with bs weights
        bprm = np.log10(self.udel*scale*bwts)
        # generate cs params with cs weights and alpha max=1
        cmax = np.minimum(self.alpha*scale, 1.0)
        cprm = np.log10(cmax*cwts)
        # join bs, cs & tau for 25 sets of random params
        rprm = np.column_stack((bprm, cprm, tprm))

        return rprm

    def random_weights(self, resample=False):
        """25 rows of random weights for bs+cs parameters"""
        # if there are no saved weights or resample is flagged, generate array
        if (not os.path.isfile(self.frnd)) or resample:
            alpha = np.ones(self.nprm)
            dims_1 = 25, self.nano
            dims_2 = 25, self.len_prm
            rwts = np.random.dirichlet(alpha, dims_1).reshape(dims_2)
            np.save(self.frnd, rwts)
        # otherwise use the saved weights
        else:
            rwts = np.load(self.frnd)
        # check dirichlet sums correctly for number of separate annos
        assert all(np.isclose(np.sum(rwts, axis=1), self.nano))

        return rwts

    @property
    def bi(self):
        """return indices of param vector for """
        return self.nba * self.nprm

    @property
    def ci(self):
        """return indices of param vector for """
        return self.bi + self.nca*self.nprm

    @property
    def nano(self):
        """return total number of annotations"""
        return self.nba + self.nca

    @property
    def len_prm(self):
        """return total length of BS+CS param vector"""
        return self.nano * self.nprm

    @property
    def npinv(self):
        """return inverse of param counts per anno"""
        return 1.0 / self.nprm

    @property
    def frnd(self):
        """return random weights file name for current configuration"""
        fstr = '{}/result/init_files/rnd_wts.bs{}.cs{}.prm{}.npy'
        return fstr.format(root_dir, self.nba, self.nca, self.nprm)


def uniform_params(umean, factor, itau, rndm=False):
    """generate randomized params summing to umean"""
    if rndm:
        rand_w = np.random.dirichlet([1]*6)
        prm = [np.log10(umean*factor*r) for r in rand_w] + [itau]
    else:
        prm = [np.log10(umean*factor / 6.0)]*6 + [itau]

    return np.array(prm)


def true_params(dist='prm'):
    """set token label for given conditions"""
    true_tau = old_cst.params[-1]
    if dist == 'prm':
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        true_prm = np.array([np.log10(mean_u / 6)] * 6 + [true_tau])
    elif dist == 'left':
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        true_prm = [np.log10(mean_u / 3)] * 3 + [-40] * 3 + [true_tau]
        true_prm = np.array(true_prm)
    elif dist == 'right':
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        true_prm = [-40] * 3 + [np.log10(mean_u / 3)] * 3 + [true_tau]
        true_prm = np.array(true_prm)
    elif dist == 'alt_1':
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        px = np.log10(mean_u / 3)
        true_prm = [px, -12.0, px, -12.0, px, -12.0] + [true_tau]
        true_prm = np.array(true_prm)
    elif dist == 'alt_2':
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        px = np.log10(mean_u / 3)
        true_prm = [-12.0, px, -12.0, px, -12.0, px] + [true_tau]
        true_prm = np.array(true_prm)
    elif 'rnd' in dist:
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        idx = int(dist.split('_')[1])
        true_prm = random_inits(idx)
    elif 'bs_cs' in dist:
        newtkn = 'ufm{}.ufmnu{}'.format(dist, prc)
        idx = int(dist.split('_')[-1])
        true_prm = random_inits(idx, 6, 6)
    else:
        newtkn = None
        true_prm = None

    return newtkn, true_prm


def random_inits(idx, nbs=6, ncs=0):
    """get random initial params from saved set"""
    true_tau = old_cst.params[-1]
    rslt = old_cst.rslt
    # generate random inits ONCE and save to file so all runs use the values
    # f_rndinit = '{}/init_files/dirich.random.inits.npy'.format(rslt)
    f_str = '{}/init_files/dirichlet.rnd.bs{}.cs{}.true.npy'
    f_rndinit = f_str.format(rslt, nbs, ncs)

    # open the saved random inits if they've already been created
    if os.path.isfile(f_rndinit):
        rndprm = np.load(f_rndinit)

    else:
        # generate random weights for bs params that sum to mean_u
        bs_wts = np.random.dirichlet([1]*nbs, 10)
        bs_prm = np.log10(mean_u * bs_wts)

        # generate random weights for cs params that sum to 0.5 (1/2 apadative)
        cs_wts = np.random.dirichlet([1]*ncs, 10)
        cs_prm = np.log10(0.5 * cs_wts)

        # 10 random sets of weights with the same tau for each set of params
        # rand_w = np.random.dirichlet([1]*6, 10)
        # jprm = np.log10(mean_u * rand_w)

        # put together an array of 10 random "true" params
        rndprm = np.column_stack((bs_prm, cs_prm, [true_tau]*10))

        # save random param distributions to file
        np.save(f_rndinit, rndprm)

    # return the specified set of random params
    return rndprm[idx]


def initial_params(idx, nbs=6, ncs=0):
    """fixed collection of initial parameter values"""
    alpha = 0.5

    # generate random inits ONCE and save to file so all runs use the values
    f_str = '{}/init_files/dirichlet.rnd.bs{}.cs{}.initial.npy'
    f_rndprm = f_str.format(old_cst.rslt, nbs, ncs)
    # f_rndprm = '{}/init_files/dirichlet.random.inits.npy'.format(old_cst.rslt)

    # open the saved random inits if they've already been created
    if os.path.isfile(f_rndprm):
        iprm = np.load(f_rndprm)

    else:
        # set initial params to single weight from 0.05, 0.5, 5 times mean_u
        ptprm = []
        for fct in (0.05, 0.5, 5):
            bprm = np.log10(mean_u * fct / nbs)
            bs_prm = [bprm] * nbs
            cprm = np.minimum(np.log10(alpha * fct / ncs), 0)
            cs_prm = [cprm] * ncs
            jprm = bs_prm + cs_prm + [init_tau]
            # jprm = [np.log10(mean_u * fct / 6.0)] * 6 + [init_tau]
            # jprm = uniform_params(mean_u, fct, init_tau)
            ptprm.append(jprm)
            for i in xrange(6):
                bprm = np.log10(mean_u*fct)
                cprm = np.minimum(np.log10(alpha*fct), 0)
                # jprm = [np.log10(prd) if j == i else -12.0 for j in xrange(6)]
                bs_prm = [bprm if j == i else -12.0 for j in xrange(nbs)]
                cs_prm = [cprm if j == i else -12.0 for j in xrange(ncs)]
                jprm = bs_prm + cs_prm + [init_tau]
                ptprm.append(jprm)

        # combine point params into array
        ptprm = np.array(ptprm)

        # get dirichlet dist of params that sum to 0.05x, 0.5x, & 5x mean_u
        rndprm = []
        for fct in (0.05, 0.5, 5):
            # 5 random sets of weights for a given factor
            # rand_w = np.random.dirichlet([1]*6, 5)
            # jprm = np.log10(mean_u * fct * rand_w)
            # rndprm.append(jprm)
            # generate random weights for bs params that sum to mean_u
            bs_wts = np.random.dirichlet([1] * nbs, 5)
            bs_prm = np.log10(fct * mean_u * bs_wts)

            # random weights for cs params that sum to 0.5 (1/2 apadative)
            cs_wts = np.random.dirichlet([1] * ncs, 5)
            cs_prm = np.minimum(np.log10(fct * alpha * cs_wts), 0)

            # put together an array of 10 random "true" params
            jprm = np.column_stack((bs_prm, cs_prm, [init_tau] * 5))
            rndprm.append(jprm)

        # save random param distributions to file
        # rndprm = np.column_stack((np.concatenate(rndprm), [init_tau]*15))
        rndprm = np.concatenate(rndprm)

        # combine point and random init params into single array and save
        iprm = np.concatenate((ptprm, rndprm))
        np.save(f_rndprm, iprm)

    # iprm = np.concatenate((ptprm, rndprm))

    return iprm[idx]




"""
LABELS AND FLAGS
"""
# set new bkgd map in cst even though it isn't needed for simulation
tkn = 'pr95.cleanrun'
b_anno = 'primate_cons95_Segments'
c_anno = 'nonsyn'
cdir = '{}_fixed_cmap'.format(c_anno)
bdir = 'std_split_pts'
prc = ''

# bdr = 'std_split_dbl_pts'
# tkn = 'dbl.precision'
# prc = '.dbl'

# bdr = 'std_split_half_pts'
# tkn = 'half.precision'
# prc = '.half'

# downsampling mask parameters
mask_id = '0001'
downsample = False


"""
SIMULATOR PARAMETERS
"""

# # double udel for each tval and simulate data from prediction
# cst.params[cst.fixed.bi:cst.fixed.bj] += np.log10(2.0)

# # halve udel for each tval and simulate data from prediction
# cst.params[cst.fixed.bi:cst.fixed.bj] -= np.log10(2.0)


"""
BOUNDARIES AND CONSTRAINTS
"""
# set bounds and inits for optimization
min_bs = None
# min_bs = np.log(0.01)

min_red = None
# min_red = 1e-02

# pi = initial_params(0, 6, 6)
# pr = random_inits(0, 6)
