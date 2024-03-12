# from scipy.optimize import LinearConstraint
from classes.dictlike import DictLike
import numpy as np

__author__ = 'davidmurphy'


class OptimizerParams(DictLike):
    """
    A class that holds the running params for and optimizer and also collects
    the results
    """
    # option settings for the different methods
    preset_options = {'COBYLA': {'disp': True, 'rhobeg': 0.1, 'tol': 1e-3,
                                 'maxiter': 10000},
                      'Powell': {'disp': True, 'ftol': 1e-05, 'maxfev': 15000,
                                 'xtol': 1e-07, 'maxiter': 4000},
                      'Nelder-Mead': {'disp': True, 'maxiter': 12000,
                                      'maxfev': 20000, 'ftol': 1e-3,
                                      'xtol': 1e-5},
                      'CG': {'disp': True, 'maxiter': 10000, 'gtol': 1e-5,
                             'norm': np.inf},
                      'BFGS': {'disp': True, 'maxiter': 10000, 'gtol': 1e-5,
                               'norm': np.inf, 'eps': None},
                      'Newton-CG': {'disp': True, 'maxiter': 10000,
                                    'xtol': 1e-5},
                      'trust-constr': {'disp': True, 'verbose':3,
                                       'maxiter': 10000,
                                       'xtol': 1e-4, 'gtol': 1e-4},
                      'TNC': {'disp': True,
                              'maxiter': 10000, 'xtol': 1e-5, 'gtol': 1e-5},
                      None: []}

    @staticmethod
    def _fcons(x, idx=0, i=0, j=0, umax=1e-6):
        """
        a constraint function for COBYLA  and SLSQP methods to use to keeping
        params in bounds on the optimization
        :param x: param vector
        :param idx: key --> (0, 1, 2) = (bs, cs, tau) constraints, respectively
        :param i: start index in param vector
        :param j: end index in param vector
        :param umax: the maximum udel for a given bs annotation
        :return: difference between max permitted values and the current values
        """
        if idx not in (0, 1, 2):
            err = 'ERROR: CONSTRAINT FUNCTION GIVEN UNKNOWN INDEX {}'
            raise IndexError(err.format(idx))
        # BACKGROUND SELECTION CONSTRAINT
        if idx == 0:
            assert (i or j)
            # max permitted udel per annotation is fixed by u_max
            return umax - np.sum(np.power(10, x[i:j]))
        # CLASSIC SWEEPS CONSTRAINT
        if idx == 1:
            assert (i or j)
            # at max 100% of alleles beneficial, i.e. alpha = 1 for each s
            return float(j - i) - np.sum(np.power(10, x[i:j]))
        # TAU CONSTRAINT
        if idx == 2:
            # tau must be within 0.1x and 10x of tau initial
            assert not (i or j)
            return (x[-1] - 0.1 * x[-1]) * (10 * x[-1] - x[-1])

    def __init__(self, method=None, bs_params=(), cs_params=(),
                 umax=7.4e-8, **kwargs):
        """
        Define a new set of optimization params for chosen method.

        :param method: a scipy.optimize.minimize method
        :param bs_params: a tuple of t-vector lengths for each bs annotation
        :param cs_params: a tuple of s-vector lengths for each cs annotation
        :param umax: maximum deleterious mutation rate allowed for bs
        """
        super(OptimizerParams, self).__init__()
        self.method = method
        self.bs_params = bs_params
        self.cs_params = cs_params
        self.umax = umax

        # "minimize" function inputs
        self.options = self.preset_options[method]
        self.target_func = None
        self.x0 = None
        self.args = None
        self.jac = None
        self.hess = None
        self.bounds = None

        # inputs for use with basinhopping algorithm
        self.niter = 100
        self.T = 1.0
        self.stepsize = 5.0

        # optimization run results
        self.runtime = 0.0
        self.fun = None
        self.status = None
        self.success = None
        self.message = None
        self.nfev = None
        self.njev = None
        self.nit = None
        self.x = None

        # initialize from keyword args
        self.setitems(items=kwargs.items())

    def update_options(self, new_options, safe=True):
        """
        Update items in option dict with "new_options" dict.

        :param new_options: options dict for current optimization method
        :param safe: flag to only set values for keys in the preset options
        :type new_options: dict
        :type safe: bool
        """
        for k, v in new_options.items():
            if k in self.options:
                self.options[k] = v
            elif (k not in self.options) and (not safe):
                self.options[k] = v
            else:
                pass

    def optimizer_results(self, result):
        """
        Write optimization results to member variables.

        :param result: an OptimizeResult object
        :type result: OptimizeResult
        """
        for k in 'fun status success message nfev njev nit x'.split():
            try:
                self[k] = result[k]
            except KeyError:
                continue

    def _bounds(self):
        """bounds on parameters being optimized"""
        bounds = []
        min_bs, max_bs = -12.0, 6.0
        if self.method == 'trust-constr':
            # for x_bs in
            pass

    def _constraints(self):
        """
        Param constraints for different optimization methods.

        List of constraint param dicts for tau, and each bs, cs annotation.
        the inequality calculated by the function 'fun' for the inputs 'args'
        must be satisfied for each param dict in the list during the
        optimization
        """
        if self.method in ['COBYLA', 'SLSQP']:

            # constraint dictionary template
            constraint_dict = {'type': 'ineq', 'fun': self._fcons}
            # gather the constraints for each parameter and selection model
            i = 0
            constraint_list = []

            # BS model constraints
            for vector_length in self.bs_params:
                j = i + vector_length
                constraint_dict['args'] = (0, i, j, self.umax)
                constraint_list.append(constraint_dict)
                i = j

            # CS model constraints
            for vector_length in self.cs_params:
                j = i + vector_length
                constraint_dict['args'] = (1, i, j)
                constraint_list.append(constraint_dict)
                i = j

            # TAU constraints
            constraint_dict['args'] = (2, 0, 0)
            constraint_list.append(constraint_dict)
            return constraint_list
        elif self.method == 'trust-constr':
            # constr = LinearConstraint()
            return ()
        else:
            return None

    # def _bhstep(self, prm):
    #     return prm * 2

    @property
    def constraints(self):
        """returns list of constraint dictionaries bs & cs annotations"""
        return self._constraints()

    @property
    def methodargs(self):
        """get the **kwargs dict to use with the optimizer function"""
        input_args = dict(fun=self.target_func, x0=self.x0, args=self.args,
                          jac=self.jac, hess=self.hess, method=self.method,
                          options=self.options, constraints=self.constraints)
        # TODO: TEMPORARY EDIT OF METHOD ARGS
        # input_args = dict(fun=self.target_func, x0=self.x0, args=self.args,
        #                   method=self.method, options=self.options,
        #                   constraints=self.constraints)
        # # input_args = dict(args=self.args,
        #                   jac=self.jac, method=self.method,
        #                   options=self.options, constraints=self.constraints)
        return input_args

    @property
    def bhargs(self):
        bhstep = lambda xa: xa * 2
        input_args = dict(func=self.target_func, x0=self.x0, niter=self.niter,
                          T=self.T, stepsize=self.stepsize, take_step=bhstep,
                          minimizer_kwargs=self.methodargs)
        return input_args