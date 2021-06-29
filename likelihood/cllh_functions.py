import numpy as np
from datetime import datetime
from sys import stdout, stderr
from parallelize import run_threads
from classes.runstruct import ChromStruct


__author__ = 'davidmurphy'


def parallel_cllh(params, cst):
    """
    Composite logLH of poly data given params for BS & CS models.

    :param params: a variable copy of free params attribute that is
                   manipulated by the optimizer and used to re-set the
                   internal attribute param values
    :param cst: chrom struct use to direct LH calculations
    :return: the negative log probability; (optional) 1st derivatives in
             param space
    :type params: numpy.ndarray
    :type cst: ChromStruct
    :rtype: float
    """

    # increment the function call counter
    cst.stat.function_calls += 1

    # # DEBUG
    # st = datetime.now()
    # message = 'starting {} thread LLH calc | '.format(cst.vars.num_cores)
    # stderr.write(message)
    # stdout.flush()

    # calculate partial sum(logLH) for each thread
    loglh_sums = run_threads(params)

    # sum the piecemeal (-1e5)-scaled logLH sums from multi-thread
    nlp = sum(loglh_sums) / cst.stat.pairs

    # # DEBUG
    # en = datetime.now()
    # message = '| LLH={}, time={}\n'.format(nlp, en - st)
    # stderr.write(message)
    # stdout.flush()

    # update internal params
    cst.params = params

    # record params, -sum(logCLH) and udel for each bs anno
    # cst.cur_params(nlp)
    # current_params = cst.cur_params(nlp)
    # stderr.write(current_params)
    # stdout.flush()

    return nlp


def serial_cllh(params, args):
    """Composite Log likelihood of neutral data given predictive LS_maps"""
    # increment the function call counter
    args.cst.stat.function_calls += 1

    # calculate partial sum(logLH) for each thread
    loglh_sum = loglh(params, *args)

    # sum the piecemeal (-1e5)-scaled logLH sums from multi-thread
    nlp = loglh_sum / args.cst.stat.pairs

    # update internal params
    args.cst.params = params

    # record params, -sum(logCLH) and udel for each bs anno
#    args.cst.cur_params(nlp)
#     current_params = args.cst.cur_params(nlp)
#     stderr.write(current_params)
#     stdout.flush()

    return nlp


def loglh(prm, cst, hom, het, u, bs, cs, s):
    """Calculate log CLH with direct inputs"""
    # predict pi across neutral sites using BS & CS maps
    pii = predicted_pi(prm, cst, u, bs, cs)

    # standard logLH calculation for the segment of neutral sites
    if s is None:
        log_lh = hom * np.log(1 - pii) + het * np.log(pii)
    else:
        log_lh = s * hom * np.log(1 - pii) + s * het * np.log(pii)

    # calculate -1e5-scaled sum inside loglh to reduce output side
    sum_llh = -1e5 * np.sum(log_lh)

    # DEBUG
    message = '. '
    stderr.write(message)
    stdout.flush()

    return sum_llh


def predicted_pi(params, cst, nu, bs, cs):
    """just return the first part of pi, dpi result from dpredicted_pi"""
    return dpredicted_pi(params, cst, nu, bs, cs)[0]
# def predicted_pi(params, cst, nu, bs, cs):
#     # # boundary conditions:
#     min_red = cst.fixed.min_red  # 1e-5
#     # min_pi0 = cst.fixed.min_pi0  # 1e-5
#     # max_pii = cst.fixed.max_pii  # 0.99
#
#     # reduction in pi from bs -> exp(dot(bs, bs_weights))
#     if bs is not None:
#         # convert bs params to udel/site for each bmap
#         uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
#         bwt = uvec / cst.fixed.u_fix
#         bsx = np.exp(np.dot(bs, bwt))
#     else:
#         bsx = np.ones(shape=len(nu))
#
#     # reduction in pi from cs -> dot(cs, cs_params)
#     if cs is not None:
#         cwt = np.power(10, params[cst.fixed.ci:cst.fixed.cj])
#         csx = np.dot(cs, cwt)
#     else:
#         csx = np.zeros(shape=len(nu))
#
#     # theta is initially mean pi scaled by mutation rate variation estimate
#     theta = nu / params[cst.fixed.ti]
#     # use meandiv instead of mean(u) -> more precise
#     # FIXED: mean_theta is not used
#     # mean_theta = cst.stat.meandiv / params[cst.fixed.ti]
#     # free param pi0 = expected het w/o selection
#     # NOTE: Remove pi0 boundary condition from original LH version
#     # pi0 = np.maximum(theta / (1.0 + theta), min_pi0)
#     pi0 = theta / (1.0 + theta)
#
#     # red = overall reduction in pi from bs & cs:
#     # FIXED: replace mean_theta with theta
#     # red = (1 + mean_theta) / (1 / bsx + csx + mean_theta)
#     red = (1 + theta) / (1/bsx + csx + theta)
#     # pii = pi after reduction from BS & CS
#     # NOTE: Remove red and pii boundary conditions from original LH version
#     # pii = np.minimum(pi0 * np.maximum(red, min_red), max_pii)
#     if min_red is not None:
#         pii = pi0 * np.maximum(red, min_red)
#     else:
#         pii = pi0 * red
#
#     return pii


def dparallel_cllh(params, cst):
    """
    Composite logLH of poly data given params for BS & CS models.

    :param params: a variable copy of free params attribute that is
                   manipulated by the optimizer and used to re-set the
                   internal attribute param values
    :param cst: chrom struct use to direct LH calculations
    :return: the negative log probability; (optional) 1st derivatives in
             param space
    :type params: numpy.ndarray
    :type cst: ChromStruct
    :rtype: float
    """

    # increment the function call counter
    cst.stat.function_calls += 1

    # # DEBUG
    # st = datetime.now()
    # message = 'starting {} thread LLH calc | '.format(cst.vars.num_cores)
    # stderr.write(message)
    # stdout.flush()

    # calculate partial sum(logLH) for each thread
    loglh_sums = run_threads(params)

    # sum the piecemeal (-1e5)-scaled logLH sums from multi-thread
    nlp = sum(ls[0] for ls in loglh_sums) / cst.stat.pairs

    # sum the parallelized gradient components if applicable
    if all(ls[1] is not None for ls in loglh_sums):
        assert cst.vars.calc_deriv
        dnlp = np.zeros(len(loglh_sums[0][1]))
        for ls in loglh_sums:
            dnlp += (np.array(ls[1]) / cst.stat.pairs)
    else:
        dnlp = None

    # # DEBUG
    # en = datetime.now()
    # message = '| LLH={}, time={}\n'.format(nlp, en - st)
    # stderr.write(message)
    # stdout.flush()

    # update internal params
    cst.params = params

    # record params, -sum(logCLH) and udel for each bs anno
    # cst.cur_params(nlp)
    current_params = cst.cur_params(nlp)
    stderr.write(current_params)
    stdout.flush()

    if cst.vars.calc_deriv:
        return nlp, dnlp
    else:
        return nlp


def dserial_cllh(params, args):
    """Composite Log likelihood of neutral data given predictive LS_maps"""
    # increment the function call counter
    args.cst.stat.function_calls += 1

    # calculate partial sum(logLH) for each thread
    loglh_sum, dloglh_sums = dloglh(params, *args)

    # sum the piecemeal (-1e5)-scaled logLH sums from multi-thread
    nlp = loglh_sum / args.cst.stat.pairs

    # gradient of function
    dnlp = np.array([dl / args.cst.stat.pairs for dl in dloglh_sums])

    # update internal params
    args.cst.params = params

    # # record params, -sum(logCLH) and udel for each bs anno
    # args.cst.cur_params(nlp)
    # current_params = args.cst.cur_params(nlp)
    # stderr.write(current_params)
    # stdout.flush()

    return nlp, dnlp


def dloglh(prm, cst, hom, het, u, bs, cs, s):
    """Calculate log CLH with direct inputs"""
    # predict pi across neutral sites using BS & CS maps
    pii, dpii = dpredicted_pi(prm, cst, u, bs, cs)

    # weighted logLH calculation with "threshold" c
    if cst.fixed.cth is not None:
        c = cst.fixed.cth  # TODO: get this from params
        # c = prm[-2]
        # get average pii vector scaled to neutral mu
        pii_avg = u * (cst.stat.meanpi / cst.stat.meandiv)
        # calculate pii, 1-pii terms including threshold components
        pii_het = ((1.0 - c) * pii) + (c * pii_avg)
        pii_hom = ((1.0 - c) * (1.0 - pii)) + (c * (1 - pii_avg))
        # calculate log likelihood with the modified terms
        log_lh = hom * np.log(pii_hom) + het * np.log(pii_het)
    else:
        log_lh = hom * np.log(1 - pii) + het * np.log(pii)

    # # standard logLH calculation for the segment of neutral sites
    # if s is None:
    #     log_lh = hom * np.log(1 - pii) + het * np.log(pii)
    # else:
    #     log_lh = s * hom * np.log(1 - pii) + s * het * np.log(pii)

    # calculate -1e5-scaled sum inside loglh to reduce output side
    sum_llh = -1e5 * np.sum(log_lh)

    if dpii is not None:
        assert cst.vars.calc_deriv
        dlog_lh = []
        for dp in dpii:
            d = hom * -dp / (1-pii) + het * dp / pii
            dlog_lh.append(d)
        # gradient of the CLH function
        dsum_llh = [-1e5 * np.sum(dl) for dl in dlog_lh]
    else:
        dsum_llh = None

    # # DEBUG
    # message = '. '
    # stderr.write(message)
    # stdout.flush()

    return sum_llh, dsum_llh


def dpredicted_pi(params, cst, nu, bs, cs):
    # # boundary conditions:
    min_red = cst.fixed.min_red  # 1e-5
    # min_pi0 = cst.fixed.min_pi0  # 1e-5
    # max_pii = cst.fixed.max_pii  # 0.99
    min_bsx = cst.fixed.min_bsx

    # reduction in pi from bs -> exp(dot(bs, bs_weights))
    if bs is not None:
        # convert bs params to udel/site for each bmap
        uvec = np.power(10, params[cst.fixed.bi:cst.fixed.bj])
        # uvec = np.array([0.0 if u < cst.fixed.u_min else u for u in uvec])
        bwt = uvec / cst.fixed.u_fix
        bsx = np.exp(np.dot(bs, bwt))
        # apply a cutoff to the composite bmap
        if min_bsx:
            bsx = np.maximum(bsx, min_bsx)
    else:
        bwt = np.zeros(shape=cst.bsparams)
        bsx = np.ones(shape=len(nu))

    # for easier use in derivative
    bsxinv = 1/bsx

    # reduction in pi from cs -> dot(cs, cs_params)
    if cs is not None:
        cwt = np.power(10, params[cst.fixed.ci:cst.fixed.cj])
        csx = np.dot(cs, cwt)
    else:
        csx = np.zeros(shape=len(nu))

    # theta is initially mean pi scaled by mutation rate variation estimate
    tau = params[cst.fixed.ti]
    theta = nu / tau

    # perform this calculation once to avoid redundancy later on
    denom = bsxinv + csx + theta

    # free param pi0 = expected het w/o selection
    pi0 = theta / (1.0 + theta)

    # TODO: can't some terms be cancelled out in advance here?
    # overall reduction in pi from bs & cs:
    red = (1 + theta) / denom

    # pi after reduction from BS & CS
    if min_red is not None:
        pii = pi0 * np.maximum(red, min_red)
    else:
        pii = pi0 * red

    # calculate the gradient if cst is flagged for derivative method
    if cst.vars.calc_deriv:
        # pre-calculate constant elements of the derivative for efficiency
        dconst = pii * (bsxinv * np.log(10)) / denom

        # partial derivative for each weight
        dpii = []
        for i in xrange(cst.bsparams):
            d = bs[:, i] * bwt[i] * dconst
            dpii.append(d)

        # partial derivative for tau
        d = -nu * (bsxinv + csx) / (tau*bsxinv + tau*csx + nu)**2.0
        dpii.append(d)
    else:
        dpii = None

    return pii, dpii
