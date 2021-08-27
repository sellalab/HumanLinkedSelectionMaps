__author__ = 'davidmurphy'


import re
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from classes.phylotree import parse_tree, load_tree, parse_exptotsub
from classes.runstruct import ChromStruct, human_autosomes, chromosome_length


class PhyloTree:
    """simplified phylogenetic tree class with nodes and branch lengths"""
    pass


def mean_tree(mod_fmt):
    """return mean phylogenetic tree from .mod format files across chroms"""
    # get the total autosomal length for weighting average branch lengths
    aut_len = 1.0 * sum(chromosome_length(c) for c in human_autosomes)

    # gather trees for each chrom using the mod file format
    meantree = Counter()
    for ch in human_autosomes:
        # get tree dict for current file
        tree = parse_tree(load_tree(mod_fmt.format(ch)))
        # weight branch lengths by chrom lengths
        weight = float(chromosome_length(ch)) / aut_len
        for (k, v) in tree.items():
            meantree[k] += float(v) * weight

    return meantree


def mean_matrix(xts_fmt):
    """return mean expected substitution matrix from .exptotsub files"""
    # initialize the dict of tree matrices with the first chromosome
    mean_mat = parse_exptotsub(xts_fmt.format(human_autosomes[0]))
    # add matrices for remaining chroms
    for ch in human_autosomes[1:]:
        # get tree dict for current file
        mat = parse_exptotsub(xts_fmt.format(ch))
        for (k, v) in mat.items():
            mean_mat[k] += v

    return mean_mat


def add_branches(add_list, out_list, mean_mat, mat):
    for br in out_list:
        if br not in add_list:
            add_list.append(br)
            mat += mean_mat[br]


def sum_branches(br_list, mean_mat):
    """sum matrices across branches in branch list"""
    # start with human branch
    mat = mean_mat['hg19']

    # make a list of outer branches that have been added
    added = []

    # create lists of outer branch configurations
    out_1 = ['hg19-panTro4']
    out_2 = out_1 + ['hg19-gorGor3']
    out_3 = out_2 + ['hg19-ponAbe2']
    out_4 = out_3 + ['hg19-nomLeu3', 'rheMac3-chlSab1']
    out_5 = out_4 + ['rheMac3-papHam1']
    out_6 = out_5 + ['rheMac3-macFas5']

    # add all relevant outer branches
    if 'gorGor3' in br_list:
        add_branches(added, out_1, mean_mat, mat)
    if 'ponAbe2' in br_list:
        add_branches(added, out_2, mean_mat, mat)
    if 'nomLeu3' in br_list:
        add_branches(added, out_3, mean_mat, mat)
    if 'chlSab1' in br_list:
        add_branches(added, out_4, mean_mat, mat)
    if 'papHam1' in br_list:
        add_branches(added, out_5, mean_mat, mat)
    if ('macFas5' in br_list) or ('rheMac3' in br_list):
        add_branches(added, out_6, mean_mat, mat)

    # add each inner branch
    for br in br_list:
        mat += mean_mat[br]

    # calculate number of matrices to divide out
    nmats = 1 + len(br_list) + len(added)

    # calculate the total number of bases across branches
    bases = mat.sum()
    # calculate total substitutions across branches
    subs = bases - mat.diagonal().sum()
    # calcuate CpG substitutions across branches
    cpg = mat[6, (4, 14)].sum()
    # adjust bases by number of trees
    bases /= float(nmats)

    # count total bases, subs and CpG adjusted subs
    b_cnt = int(bases)
    s_cnt = int(0.5 * subs)
    c_cnt = int(0.5 * (subs - cpg))

    # return b_cnt, s_cnt, c_cnt
    return 1.0 * c_cnt / b_cnt, 1.0 * s_cnt / b_cnt


def main():
    fmt = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast/' \
          'primate_neutral/{}.neutral.exptotsub'
    mm = mean_matrix(fmt)
    # blist = ['hg19', 'panTro4', 'rheMac3', 'ponAbe2', 'macFas5', 'nomLeu3']
    # dist = sum_branches(blist, mm)
    astats = '/Users/davidmurphy/GoogleDrive/linked_selection/data/phast/' \
             'align_stats_combined.txt'
    lpct, lcdist, lsdist, lnum = [], [], [], []
    with open(astats, 'r') as f:
        for line in f:
            sp, _, pct = line.split()
            slist = sp.split('_')[1:]
            cdist, sdist = sum_branches(slist, mm)
            lpct.append(float(pct))
            lnum.append(len(slist))
            lcdist.append(cdist)
            lsdist.append(sdist)
            print '{} {} {} {} {}'.format(sp, pct, cdist, sdist, cdist*float(pct))

    lpct, lcdist, lsdist, lnum = [np.array(a) for a in lpct, lcdist, lsdist, lnum]
    # plt.scatter(ldist, lpct)
    # plt.xlabel('total branch length')
    # plt.ylabel('percent neutral alignment')
    # plt.show()
    # plt.scatter(lnum, ldist)
    # plt.xlabel('aligned species')
    # plt.ylabel('total branch length')
    # plt.show()
    # plt.scatter(lcdist, lpct, color='darkorange', alpha=0.8, label='nonCpG')
    plt.scatter(lsdist, lpct, color='purple', alpha=0.8, label='all')
    # for xmx in np.arange(0.018, 0.16, 0.01):
    #     xi = np.arange(0.001, 1.5*xmx, 0.001)
    #     yi = xmx / xi
    #     plt.plot(xi ,yi, color='k', alpha=0.75)
    plt.xlabel('total branch length')
    plt.ylabel('percent neutral alignment')
    plt.ylim(0.8, 1)
    plt.xlim(0.015, 0.17)
    plt.legend()
    fsave = '/Users/davidmurphy/GoogleDrive/linked_selection/result/' \
            'final_files/branch_lengths.png'
    plt.savefig(fsave, dpi=256)
    plt.close()
    # plt.show()


if __name__ == '__main__':
    main()