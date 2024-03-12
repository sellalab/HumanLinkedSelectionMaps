__author__ = 'davidmurphy'


import os
import gzip
import pickle
import sys


final_dir = r'/Users/MURPHYD/Dropbox (OMRF)/final_files'


# import seaborn
# import numpy as np
# import matplotlib.pyplot as plt
# from data_processing.functions import count_all_pairs, calc_pi
# from classes.runstruct import ChromStruct, root_dir
# from figures.common_functions import format_panels, cst_from_fldr

# average over 1,2,5,10mb with all neutral sites
# mark chrom ends, do manhattan style plot

# #%%
# fpickle = '/Users/MURPHYD/Desktop/variants.pickle.gz'
# ch22 = pickle.load(gzip.open(fpickle, 'r'))
# positions = sorted(ch22.keys())
# alt_counts = np.array([ch22[p][1] if 1 in ch22[p] else 0 for p in positions])
# pi_vals = calc_pi(216, alt_counts)
#
# #%%
# px, pii = [], []
# dist = int(len(positions)/500)
# for i in range(0, len(positions), dist):
#     px.append(np.mean(positions[i:i+dist]))
#     pii.append(np.mean(pi_vals[i:i+dist]))
# plt.figure()
# plt.scatter(px, pii)
# plt.show()


def build_mock_frq_file_from_neutral_coal_simulations(chrom):
    """
    build a mock .frq file using the variation data from Will's
    simulations of the neutral coalesent
    """
    # load the pickle data
    coaldir = final_dir + '/chromosomeARGS'
    coalpickle = coaldir + '/tree_chromosome_{}/variants.pickle.gz'.format(chrom)
    chdict = pickle.load(gzip.open(coalpickle, 'r'))
    pos = sorted(chdict.keys())  # get positions in order

    # iterate through all of the position and create mock .frq file lines
    # mockdir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/snps/frq/mockYRI'
    mockdir = final_dir + '/mockYRI'
    os.makedirs(mockdir, exist_ok=True)
    mockfile = mockdir + '/chr{}.mockYRI.all.phase3.frq.count.gz'.format(chrom)
    linefmt = '{}\t{}\t2\t216\t{}\t{}\n'  # format for .frq lines
    # 11      73015   2       216     C:216   T:0
    with gzip.open(mockfile, 'w') as f:
        for p in pos:
            # note: there were a small count of site with only ref
            if 1 not in chdict[p]:
                continue
            # for each simulated position, use A/T (doesn't matter) for alleles
            altcount = 'A:{}'.format(chdict[p][1])
            refcount = 'T:{}'.format(216 - chdict[p][1])
            line = linefmt.format(chrom, int(p+1), refcount, altcount)
            f.write(line.encode('utf-8'))


def main():
    for c in range(1, 23):
        build_mock_frq_file_from_neutral_coal_simulations(c)
#     if len(sys.argv) != 2:
#         print('usage: neutral_coalescent <chrom>')
#         exit(1)
#     chrom = sys.argv[1]
#     build_mock_frq_file_from_neutral_coal_simulations(chrom)
main()

if __name__ == '__main__':
    main()
