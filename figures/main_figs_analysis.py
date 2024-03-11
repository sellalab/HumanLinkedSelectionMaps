__author__ = 'davidmurphy'


import os
import subprocess
import numpy as np
from sys import argv
from classes.runstruct import ChromStruct, root_dir
from data_processing.functions import harmonic_mean


def main():
    # folder name for the main CADD 6% results:
    fldr = 'cadd94_gmask_v1.6_without_bstat'

    # use subprocess to call each analysis script:

    # run analysis for fig 2A observed vs. predicted chromosome 1
    fig2a_call = '{}/lsm_python/likelihood/chrom_map.py {}'.format(root_dir, fldr)
    subprocess.call(fig2a_call, shell=True)
    # run calculate R^2 analysis for main fig 2B
    fig2b_call = '{}/lsm_python/likelihood/calc_rsq.py {}'.format(root_dir, fldr)
    subprocess.call(fig2b_call, shell=True)
    # run collate diversity analysis around nonsynonymous substitutions (fig 3)
    fig3_call = '{}/lsm_python/results/collate_diversity.py {} 5e-3 nonsyn'.format(root_dir, fldr)
    subprocess.call(fig3_call, shell=True)
    # run data sorted by B value anaylsis (fig 5)
    fig5_call = '{}/lsm_python/likelihood/sort_pred.py {}'.format(root_dir, fldr)
    subprocess.call(fig5_call, shell=True)


if __name__ == '__main__':
    main()
