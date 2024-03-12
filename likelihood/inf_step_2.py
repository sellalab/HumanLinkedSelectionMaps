__author__ = 'davidmurphy'


import os
import numpy as np
from sys import argv
from runinf import run_inference
from classes.runstruct import ChromStruct, root_dir
from data_processing.functions import harmonic_mean


def run_mean_best(folder_name):
    fdir = root_dir + '/result/final_files/{}/'.format(folder_name)
    # collect results files and create ChromStructs
    rlst = []
    flst = []
    for f in os.listdir(fdir):
        if f.endswith('.txt'):
            fname = fdir + f
            flst.append(fname)
            rlst.append(ChromStruct(chrom='chr1', init=fname))

    # assign indices to each runstruct
    ilh = [(i, r.stat.best_lh) for (i, r) in enumerate(rlst)]
    # get indices of the top 3 best LH runs (after sorting on LH)
    top3 = [i for (i, r) in sorted(ilh, key=lambda x: x[1])[:3]]

    # # get the average params
    # avg_params = np.average([rlst[i].params for i in top3], axis=0)

    # get the harmonic mean of the params
    prm_group = np.array([rlst[i].params for i in top3])
    avg_params = []
    for i in range(len(rlst[0].params)):
        m = harmonic_mean(prm_group[:,i])
        avg_params.append(m)

    # create new ChromStruct (use any result file for init)
    cst = ChromStruct(chrom='chr1', init=flst[0])
    cst.params = avg_params
    cst.stat.init_params = cst.params

    # rename final text file
    fsave = cst.init_file.replace('.final.', '.final.composite.')
    cst.txt_file = fsave

    run_inference(cst, parallel=True)


def main():
    # folder name for the main CADD 6% results:
    fldr = 'cadd94_gmask_v1.6_without_bstat'
    run_mean_best(fldr)

    # if os.getcwd().startswith('/Users/davidmurphy/'):
    #     fldr = 'nff2'
    #     run_mean_best(fldr)
    #
    # else:
    #     if len(argv) != 2:
    #         print('usage: mean_best <folder_name>')
    #         exit(1)
    #     else:
    #         run_mean_best(argv[1])


if __name__ == '__main__':
    main()