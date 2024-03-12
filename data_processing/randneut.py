from classes.runstruct import ChromStruct
from data_tools import randint_unique
from sys import argv
import numpy as np


__author__ = 'davidmurphy'


# def rand_neut_file(fneut, n):
#     suff = '.random.sample.n_{}.npz'.format(n)
#     fsave = fneut.replace('.npz', suff).replace('/mask/', '/rand/')
#     return fsave


def select_sites(n_sites, chrom, tok, save=False):
    """select n neutral sites randomly from neutral mask file"""
    # initialize chrom struct to link to file paths
    cst = ChromStruct(chrom, tkn=tok)
    fmask = cst.neut_masks

    # get positions of all neutral sites in the neutmask file
    neut_sites = np.where(np.load(fmask)['neutmask'] > 0)[0]

    # select sites using random indices, shift coordinate system up by 1
    ridx = randint_unique(n_sites, neut_sites.size-1)
    sample = neut_sites[ridx] + 1

    # either save the file or return the random sample
    if save:
        # new_suffix = '.random.sample.n_{}.npz'.format(n_sites)
        # fsave = fmask.replace('.npz', new_suffix).replace('/mask/', '/rand/')
        fsave = cst.rand_file(n_sites)
        np.savez_compressed(fsave, randneut=sample)
    else:
        return sample


def main():
    if len(argv) != 3:
        print 'usage: random_neutral_sites <n_sites> <chrom>'
        exit(1)

    n_sites = int(eval(argv[1]))
    chrom = argv[2]
    tok = 'pr95.cleanrun'
    select_sites(n_sites, chrom, tok, True)


if __name__ == '__main__':
    main()
    # select_sites(50000, 'chr22', 'pr95.cleanrun', True)
