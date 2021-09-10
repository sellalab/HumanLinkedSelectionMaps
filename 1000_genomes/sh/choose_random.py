#!/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python
import numpy as np
from sys import argv

__author__ = 'davidmurphy'


def main():
    """Pick n random individuals from a given population"""
    if len(argv) != 4:
        print 'usage: choose_random <N> <POP> <ID_FILE>'
        exit(1)

    # root dir
    root = '/ifs/data/c2b2/gs_lab/dam2214/pyLS/snps'

    # set the sample size and population ID
    n = argv[1]
    pop_id = argv[2]

    # load the individual IDs and associated population name for each individual
    ind, pop = np.loadtxt(argv[3], usecols=(0, 1), dtype=str).T

    # get the subset of IDs for each population
    subpop = ind[pop == pop_id]
    # make a random index to select 81 of the subset of IDs
    rand_ind = np.random.choice(subpop, size=81, replace=False)

    # save the random subset in a file containing N and population ID
    np.savetxt('{}/pops/{}.random.{}.IDs.txt'.format(root, pop_id, n), rand_ind, fmt='%s')


if __name__ == '__main__':
    main()
