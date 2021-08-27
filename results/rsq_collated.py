import os
from sys import argv, stderr, stdout
from functions import rsquared_function
from results.mapfunctions import *

__author__ = 'davidmurphy'


def main_remote():
    if len(argv) == 4:
        init = argv[1]
        rst = RunStruct(init=init)
        scale = float(argv[2])
        out_file = argv[3]
        collated_file = '{}/result/collate/{}_{}_{}.npy'.format(rst.root, rst.token, rst.focal, rst.neut)
        rsq_collated(collated_file=collated_file, scale=scale, out_file=out_file)
    else:
        print 'usage: rsq_collated <init_file> <scale> <out_file>'


def main_local():
    root = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/'
    cf = root + 'result/collate/cds-mid/prosimian-cons-95-segs_cdsMid_YRI_neutral.npy'
    f = root + 'result/rsq_collated/primate_cdsMid_YRI_neutral.txt'
    for scale in [1e-4, 1e-3, 1e-2]:
        rsq_collated(collated_file=cf, scale=scale, out_file=f)


def rsq_collated(collated_file, scale, out_file):
    """
    Calculate rsquared across different window sizes for the collate diversity data
    :param collated_file: the file to perform the calculations on
    :param scale: window size for averaging data
    :param out_file: file to write [window, rsq] results to
    """
    win, pi, div, pred, cnts = np.load(collated_file).T
    # pi and pred are divided by counts to get mean value/window. then they are normalized to 1
    pi /= cnts
    pred /= cnts
    # normalize to 1
    # pi /= np.mean(pi)
    # pred /= np.mean(pred)

    # binned averages of pi, pred
    # from analyze_inf_results.moving_avg import moving_avg
    # bins, pib = moving_avg(np.column_stack([win, pi]), window=scale, overlap=0.5 * scale).T
    # bins, prb = moving_avg(np.column_stack([win, pred]), window=scale, overlap=0.5 * scale).T
    # bns, pib = moving_avg(np.column_stack([win, pi]), window=scale, overlap=scale).T
    # bns, prb = moving_avg(np.column_stack([win, pred]), window=scale, overlap=scale).T
    # bins, pib = binned_function(win, pi, scale, step=0.5 * scale)
    # bins, prb = binned_function(win, pred, scale, step=0.5 * scale)
    bins, pib = binned_function(win, pi, scale)
    bins, prb = binned_function(win, pred, scale)

    # remove nan windows prior to rsq calc
    mask = (np.isfinite(pib) & np.isfinite(prb))
    rsq = rsquared_function(xobs=pib[mask], yobs=prb[mask])
    # [window, rsq] results format is used
    resultstring = '{:.2e} {}\n'.format(scale, rsq)
    stderr.write(resultstring)
    stdout.flush()
    with open(out_file, 'a') as f:
        f.write(resultstring)

    return None


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
