__author__ = 'davidmurphy'

import numpy as np
from sys import argv
from likelihood.cllh_functions import serial_cllh
from precalc.lh_inputs import prepare_inputs
from classes.runstruct import ChromStruct, root_dir


def init_cllh(cst, params):
    cst.vars.num_cores = 1
    args = prepare_inputs(cst)
    return serial_cllh(params, args)


def main():
    fdir = '{}/result/final_files/sims'.format(root_dir)
    fn = '00-YRI.pr95.clustinf.initial.BS1.6.CS0.0.180408214759.final.txt'
    fold = '{}/{}'.format(fdir, fn)
    old_cst = ChromStruct(chrom='chr1', init=fold)
    old_cst.stat.calc_stats(old_cst)
    mean_u = old_cst.stat.utot[0]

    if root_dir.startswith('/Users/davidmurphy'):
        finit = fold
        inf_cst = old_cst

    else:
        if len(argv) != 2:
            print 'usage: print_cllh <init>'
            exit(1)
        finit = argv[1]
        inf_cst = ChromStruct('chr1', init=finit)

    # if 'ufmprm' in finit:
    #     params = [np.log10(mean_u / 6)] * 6 + [old_cst.params[-1]]
    # elif 'ufmleft' in finit:
    #     params = [np.log10(mean_u / 3)] * 3 + [-40] * 3 + [old_cst.params[-1]]
    # else:
    #     params = [-40] * 3 + [np.log10(mean_u / 3)] * 3 + [old_cst.params[-1]]
    #
    # params = np.array(params)

    print init_cllh(inf_cst, inf_cst.init_params)
    print init_cllh(inf_cst, inf_cst.stat.best_params)

    return None


if __name__ == '__main__':
    main()