import numpy as np
from sys import argv, stderr, stdout
from datetime import datetime
from multiprocessing import cpu_count
from scipy.optimize import minimize as minimize
from scipy.optimize import basinhopping as basinhopping
from precalc.lh_inputs import prepare_inputs
from classes.runstruct import ChromStruct, OptimizerParams, root_dir
from parallelize import init_threads, set_fixed_params, set_state
from data_processing.data_tools import time2str, str2time
from cllh_functions import serial_cllh, dserial_cllh, \
     dparallel_cllh, parallel_cllh, loglh, dloglh


__author__ = 'davidmurphy'


def evaluate(cst, grad, test_prm):
    # get inputs for optimization run
    inputs = prepare_inputs(cst)
    params = cst.params
    last_params = np.copy(params)
    last_nlp, last_dnlp = dserial_cllh(params, inputs)

    for prm in xrange(7):
        print 'param {}'.format(prm+1)
        # for x in np.logspace(-6, -2, 5):
        for x in [1e-6]:
            params[prm] -= params[prm] * x
            deltap = last_params[prm] - params[prm]
            nlp, dnlp = dserial_cllh(params, inputs)
            print (last_nlp - nlp) / deltap, last_dnlp[prm]
            last_nlp, last_dnlp, last_params = nlp, dnlp, np.copy(params)
        print '----\n'

    return None


def basin_hopper(cst):
    # record run start time format=YYMMDDHHMMSS
    cst.record_timestamp()
    # affix start time to init_file to
    cst.save()

    # switch off the 'complete' flag before entering optimization
    cst.vars.complete = False

    # loop through optimization methods
    for optimizer_id in cst.optimizers:
        # store start/end times for each optimization method
        start_time = datetime.now()

        # get optimizer data structure for 'optimizer_id' to parameterize
        optimizer = cst[optimizer_id]

        # set initial parameter to current cst.params
        optimizer.x0 = cst.params

        # run optimizers in sequence and use each optimizer result as
        # the set of initial params for the following optimizer
        # optimization = basinhopping(**optimizer.bhargs)
        optimization = minimize(**optimizer.methodargs)

        # process OptimizeResult object
        optimizer.optimizer_results(result=optimization)

        # get final params, use as x0 for next method
        cst.params = optimizer.x

        # update best method, likelihood and param values
        if optimizer.fun < cst.stat.best_lh:
            cst.stat.best_lh = optimizer.fun
            cst.stat.best_method = optimizer.method
            cst.stat.best_params = optimizer.x

        # save the method run time and cumulative run time
        optimizer.runtime = time2str(datetime.now() - start_time)
        hours = str2time(optimizer.runtime).total_seconds() / 3600.0
        cst.stat.total_time += round(hours, 2)

    # write the final output to the optimization log
    if len(cst.cached_params) > 0:
        cst.log_params()

    # mark the run as complete and calculate final stats
    cst.vars.complete = True
    cst.stat.calc_stats(cst)
    cst.save(txt_file=cst.final_file)


def optimize_1(cst, optimizer):
    """
       Notes on some of the optimization methods used:
       -- Nelder-Mead is fast and takes large steps but it cannot be bounded
       so it may return negative weights
       -- In practice we fix the negative weights in the log_likelihood
       function and call them either 0 or min_w
       -- Truncated Newton is close to the 'active-set' option used in MATLAB
       """
    # switch off the 'complete' flag before entering optimization
    cst.vars.complete = False

    # store start/end times for each optimization method
    start_time = datetime.now()

    # set initial parameter to current cst.params
    optimizer.x0 = cst.params

    # run optimization with the selected optimizer
    optimization = minimize(**optimizer.methodargs)

    # process OptimizeResult object
    optimizer.optimizer_results(result=optimization)

    # get final params, use as x0 for next method
    cst.params = optimizer.x

    # update best method, likelihood and param values
    if optimizer.fun < cst.stat.best_lh:
        cst.stat.best_lh = optimizer.fun
        cst.stat.best_method = optimizer.method
        cst.stat.best_params = optimizer.x

    # save the method run time and cumulative run time
    optimizer.runtime = time2str(datetime.now() - start_time)
    hours = str2time(optimizer.runtime).total_seconds() / 3600.0
    cst.stat.total_time += round(hours, 2)

    # return cst


def parallel_optimization_1(cst):
    # get number of cores
    n = cst.vars.num_cores

    # set functions for optimization and parallelized threads
    f_optimize, f_thread = dparallel_cllh, dloglh

    # start processes prior to loading data
    # init_threads(nproc=n, func=f_thread)

    # get inputs for optimization run
    inputs = prepare_inputs(cst)

    # loop through optimization methods
    for optimizer_id in cst.optimizers:

        # get optimizer data structure for 'optimizer_id' to parameterize
        optimizer = cst[optimizer_id]

        # DEBUG
        optimizer.options['maxiter'] = 10
        optimizer.options['maxfev'] = 10

        # TODO: internalize this step somewhere else at some point
        if optimizer.method in ['BFGS', 'CG']:
            jac = True
            cst.vars.calc_deriv = True
        else:
            jac = None
            cst.vars.calc_deriv = False

        # set keyword arguments for optimizer function
        cst.optimizer_inputs(func=f_optimize, args=(cst,), jac=jac)

        # set/reset the inputs for each sub-process
        init_threads(nproc=n, func=f_thread)
        set_fixed_params(inputs)

        optimize_1(cst, optimizer)

    # write the final output to the optimization log
    if len(cst.cached_params) > 0:
        cst.log_params()

    # mark the run as complete and calculate final stats
    cst.vars.complete = True
    cst.stat.calc_stats(cst)
    cst.save(txt_file=cst.final_file)


def serial_optimization(cst, grad=False):
    assert cst.vars.num_cores == 1

    # get inputs for optimization run
    inputs = prepare_inputs(cst)

    # set optimizer params and data input args
    if grad:
        cst.optimizer_inputs(func=dserial_cllh, args=(inputs,))
    else:
        cst.optimizer_inputs(func=serial_cllh, args=(inputs,))

    # run the optimization with a single processor
    optimize(cst)


def optimize(cst):
    """
    Run optimization with one or more algorithms.
    :param cst: structure organizing all inputs and outputs for optimization.
    :type cst: ChromStruct
    """
    # if using simulation data with "true" params, record their CLH
    if cst.stat.true_params is not None:
        cst.stat.true_clh = dparallel_cllh(cst.stat.true_params, cst)

    # record the likelihood of initial params
    cst.stat.init_clh = dparallel_cllh(cst.stat.init_params, cst)

    # record run start time format=YYMMDDHHMMSS
    cst.record_timestamp()
    # # save initial file
    # cst.save()

    # switch off the 'complete' flag before entering optimization
    cst.vars.complete = False

    # loop through optimization methods
    for optimizer_id in cst.optimizers:
        # store start/end times for each optimization method
        start_time = datetime.now()

        # get optimizer data structure for 'optimizer_id' to parameterize
        optimizer = cst[optimizer_id]
        assert isinstance(optimizer, OptimizerParams)

        # set initial parameter to current cst.params
        optimizer.x0 = cst.params

        # run optimizers in sequence and use each optimizer result as
        # the set of initial params for the following optimizer
        # optimization = basinhopping(**optimizer.bhargs)
        optimization = minimize(**optimizer.methodargs)

        # process OptimizeResult object
        optimizer.optimizer_results(result=optimization)

        # get final params, use as x0 for next method
        cst.params = optimizer.x

        # update best method, likelihood and param values
        if optimizer.fun < cst.stat.best_lh:
            cst.stat.best_lh = optimizer.fun
            cst.stat.best_method = optimizer.method
            cst.stat.best_params = optimizer.x

        # save the method run time and cumulative run time
        optimizer.runtime = time2str(datetime.now() - start_time)
        hours = str2time(optimizer.runtime).total_seconds() / 3600.0
        cst.stat.total_time += round(hours, 2)

    # write the final output to the optimization log
    if len(cst.cached_params) > 0:
        cst.log_params()

    # mark the run as complete and calculate final stats
    cst.vars.complete = True
    cst.stat.calc_stats(cst)
    cst.save(txt_file=cst.final_file)


def parallel_optimization(cst, use_basinhopping=False):
    # get number of cores
    n = cst.vars.num_cores

    # set functions for optimization and parallelized threads
    f_optimize, f_thread = dparallel_cllh, dloglh

    # start processes prior to loading data
    init_threads(nproc=n, func=f_thread)

    # get inputs for optimization run
    inputs = prepare_inputs(cst)

    # set the inputs for each sub-process
    set_fixed_params(inputs)

    # set keyword arguments for optimizer function
    cst.optimizer_inputs(func=f_optimize, args=(cst,))

    # if not root_dir.startswith('/Users/davidmurphy'):
    #     # run the optimization on the threads
    #     if use_basinhopping:
    #         # cst = basin_hopper(cst)
    #         optimize(cst)
    #
    #     else:
    #         # msg = 'running optimization\n'.format(str(cst.params))
    #         # stderr.write(msg)
    #         # stdout.flush()
    #         optimize(cst)

    optimize(cst)


def run_inference(cst, parallel=True):
    # run in serial on a single core
    if not parallel:
        cst.vars.num_cores = 1
        pass

    # run threaded optimization on all available cores
    else:
        cst.vars.num_cores = cpu_count()
        assert cst.vars.num_cores > 1
        parallel_optimization(cst)

    # return cst
