import os
import numpy as np
from sys import argv, stdout, stderr
from precalc.lh_inputs import load_saved, adjust_arrays
from likelihood.cllh_functions import predicted_pi
from classes.runstruct import ChromStruct, root_dir, izip
from data_processing.functions import rsquared_function, swap_root


__author__ = 'davidmurphy'


def rsq_2(cst, scale, nt, nu, bs, cs, dv, sg):
    # format arrays for the predicted_pi function
    bs, _, nu, nt, dv, _, ms = adjust_arrays(cst, bs, cs, nu, nt, dv, None)

    # def adjust_arrays(cst, bs, cs, nu, nt, dv, pl):
    # bs, cs, nu, nt, dv, pl, ms

    # get cumulative genomic positions for each upper segment edge
    cum_pos = np.cumsum(sg)
    # get windows tiling across the cumulative positions
    windows = np.arange(scale, cum_pos[-1], scale)

    # get indicies to sort segments into windows based on cumulative position
    jdx = np.searchsorted(cum_pos[ms], windows)
    idx = np.concatenate(([0], jdx[:-1]))

    # get the vector of predicted pi values per segment
    pred = predicted_pi(cst.params, cst, nu, bs, None)

    # observed, predicted, number of data points per Mb
    obs, prd, num = [], [], []
    for (i, j) in izip(idx, jdx):
        if j > i:
            # calculate summaries on sites sorted into each window
            pi = nt[i:j, 1].sum() / nt[i:j].sum()
            pr = np.average(pred[i:j], weights=dv[i:j])
            n = dv[i:j].sum()
            # append summary stats to lists
            obs.append(pi)
            prd.append(pr)
            num.append(n)

    # convert lists to arrays
    obs, prd, num = map(np.array, (obs, prd, num))
    msk = (num > 0.2*np.mean(num))  # ...
    rsq = rsquared_function(obs[msk], prd[msk])

    return rsq


def rsq_maps(mst, scale):
    """
    This function is used to calculate rsquared for observed vs. predicted
    diversity levels at a given spatial scale and for a given set of inferred
    parameters and data. The scale-rsquared values are printed to stderr and
    also sent to a file.
    :param mst: MapStruct containing params and arrays
    :param scale: scale size in bp on which to calculate rsquared
    """
    obs, pred, sites = [], [], []
    jdx = mst.window_index(window=scale)
    idx = np.concatenate([[0], jdx[:-1]])
    prediction = mst.prediction()
    n = mst.nsites
    for (i, j) in izip(idx, jdx):
        if j > i:
            pi = np.sum(mst.het[i:j]) / np.sum(mst.het[i:j]+mst.hom[i:j])
            pr = np.average(prediction[i:j], weights=n[i:j])
            st = np.sum(n[i:j])
        else:
            pi = np.nan
            pr = np.nan
            st = 0

        obs.append(pi)
        pred.append(pr)
        sites.append(st)

    # convert lists to arrays
    obs, pred, sites = map(np.array, [obs, pred, sites])
    # mask nan sites
    ni = np.isfinite(obs) & np.isfinite(pred)
    obs, pred, sites = obs[ni], pred[ni], sites[ni]

    # mask nan segments and segments where site density < 20% of mean
    msk = (sites > (np.mean(sites) * 0.2))

    # return calculated r^2 on the masked data
    return rsquared_function(xobs=obs[msk], yobs=pred[msk])


def main_remote():
    if len(argv) != 2:
        print 'usage: rsq_maps <init>'
        exit(1)

    # load completed inference results
    cst = ChromStruct('chr1', init=argv[1])
    # new, sim = map(eval, argv[2:])

    # load arrays
    #     return sg, bs, cs, nu, nt, dv, pl
    sg, bs, cs, nu, nt, dv, _ = load_saved(cst)

    # get R^2 over a set of log scaled windows
    scales = np.concatenate((np.logspace(4, 6, 20),
                             np.logspace(6.11111111, 6.555555555, 5)))

    # write header for tabular stderr printout
    stderr.write('scale R^2\n')
    stdout.flush()

    # clear any existing R^2 values in stats
    cst.stat.rsq_values = []

    # record in results in stats and print to stderr
    for scl in scales:
        rsquared = rsq_2(cst, scl, nt, nu, bs, cs, dv, sg)
        cst.stat.rsq_values.append((scl, rsquared))
        stderr.write('1e{:.3f} {:.4f}\n'.format(np.log10(scl), rsquared))
        stdout.flush()

    # re-save cst
    cst.save()


def main_local():
    fdir = '{}/result/final_files/sims'.format(root_dir)
    csfiles = ['{}/{}'.format(fdir, f) for f in os.listdir(fdir) if
               f.endswith('.txt')]
    for cf in csfiles:
        swap_root(cf)
    # load completed inference results
    cst = ChromStruct('chr1', init=csfiles[5])
    # new, sim = map(eval, argv[2:])

    # load arrays
    #     return sg, bs, cs, nu, nt, dv, pl
    sg, bs, cs, nu, nt, dv, _ = load_saved(cst)

    # # get R^2 over a set of log scaled windows
    # scales = np.concatenate((np.logspace(4, 6, 20),
    #                          np.logspace(6.11111111, 6.555555555, 5)))

    # write header for tabular stderr printout
    stderr.write('scale R^2\n')
    stdout.flush()

    # clear any existing R^2 values in stats
    cst.stat.rsq_values = []

    # record in results in stats and print to stderr
    for scl in [1e6]:
        rsquared = rsq_2(cst, scl, nt, nu, bs, cs, dv, sg)
        cst.stat.rsq_values.append((scl, rsquared))
        stderr.write('1e{:.3f} {:.4f}\n'.format(np.log10(scl), rsquared))
        stdout.flush()

    # re-save cst
    # cst.save()


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
