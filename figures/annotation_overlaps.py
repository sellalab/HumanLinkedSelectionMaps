__author__ = 'davidmurphy'

import os
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from classes.runstruct import ChromStruct, root_dir, human_autosomes
from figures.common_functions import format_panels


#%%
def overlap_two_annos(cst, anno1, anno2):
    """count overlapping and nonoverlapping sites in two annotations"""
    assert isinstance(cst, ChromStruct)

    # set two empty arrays for each annotation set
    a_1 = np.zeros(shape=cst.chlen, dtype=bool)
    a_2 = np.zeros(shape=cst.chlen, dtype=bool)

    # flip 0s to 1s in each annotation array using bed files
    f_1 = cst.bs_target(anno1)
    with open(f_1, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            start, end = map(int, line.split()[1:])
            a_1[start:end] = True

    f_2 = cst.bs_target(anno2)
    with open(f_2, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            start, end = map(int, line.split()[1:])
            a_2[start:end] = True

    # get the sum of sites for each annotation type and their overlap
    c_1 = np.sum(a_1)
    c_2 = np.sum(a_2)
    overlap = np.sum(a_1 & a_2)

    return c_1, c_2, overlap


cst = ChromStruct('chr1')
a1 = 'fish_cons94_gmask'
a2 = 'exon'
overlap_two_annos(cst, a1, a2)


#%% OVERLAP PHASTCONS AND CADD
def overlap_genic(cst, anno, compareto='genic'):
    """count overlapping and nonoverlapping in genic annos and other anno"""
    assert isinstance(cst, ChromStruct)

    # set two empty arrays for each annotation set
    a_1 = np.zeros(shape=cst.chlen, dtype=bool)

    # flip 0s to 1s in each annotation array using bed files
    f_1 = cst.bs_target(anno)
    with open(f_1, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            start, end = map(int, line.split()[1:])
            a_1[start:end] = True

    # get arrays from each annotation
    if compareto == 'genic':
        genic_annos = ['cds', 'utr', 'tss', 'peri', 'intron']
    else:
        genic_annos = ['cds', 'cCRE_PLS_filtered', 'cCRE_ELS_filtered',
                       'cCRE_CTCF_filtered', 'cCRE_H3K4me3_filtered']
    genic_arras = [np.zeros(shape=cst.chlen, dtype=bool) for _ in genic_annos]
    # make an array of ones for intergenic
    inter_arra = np.ones(shape=cst.chlen, dtype=bool)
    for (an, ar) in zip(genic_annos, genic_arras):
        f_in = cst.bs_target(an)
        with open(f_in, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                start, end = map(int, line.split()[1:])
                ar[start:end] = True
                # remove all genic annotations from intergenic
                inter_arra[start:end] = False

    # get the sum of sites for each annotation type and their overlap
    c_1 = np.sum(a_1)
    genic_counts = [np.sum(ar) for ar in genic_arras]
    # include the intergenic array
    genic_counts.append(np.sum(inter_arra))
    overlaps = [np.sum(a_1 & ar) for ar in genic_arras]
    # include the overlap to intergenic
    overlaps.append(np.sum(a_1 & inter_arra))

    return c_1, genic_counts, overlaps


def save_overlap_analysis(annotation, compareto='genic'):
    final_dir = root_dir + '/result/final_files'
    n_cons = 0
    num_annos = 6 if compareto == 'genic' else 6
    total_counts = [0] * num_annos
    overl_counts = [0] * num_annos
    for ch in human_autosomes:
        cst = ChromStruct(ch)
        c1, cg, co = overlap_genic(cst, annotation)
        n_cons += c1
        for i in xrange(num_annos):
            total_counts[i] += cg[i]
            overl_counts[i] += co[i]

    if compareto == 'genic':
        f_save = final_dir + '/sfigs/overlap_{}.txt'.format(annotation)
    else:
        f_save = final_dir + '/sfigs/ccre_overlap_{}.txt'.format(annotation)
    with open(f_save, 'w') as f:
        for (tot, ovr) in zip(total_counts, overl_counts):
            line = '{} {}\n'.format(tot, ovr)
            f.write(line)


a1 = 'fish_cons94_new'
a2 = 'cadd94_gmask_v1.6_without_bstat'
# save_overlap_analysis(a1, 'conserved')
# save_overlap_analysis(a2, 'CADD')
# save_overlap_analysis(a2, 'genic')
save_overlap_analysis(a1, 'cCRE')
save_overlap_analysis(a2, 'cCRE')


#%% PLOT OVERLAPS (VERSION USED IN S4)
def plot_combine_cadd_cons(anno_type='genic'):
    a1 = 'fish_cons94_new'
    a2 = 'cadd94_gmask_v1.6_without_bstat'
    frac_in = []
    frac_of = []

    # used for range operations on lists
    if anno_type == 'genic':
        anno_num = 6
    else:
        anno_num = 6

    # load data, calculate fractions
    for annotation in [a1, a2]:
        # set path to saved data
        final_dir = root_dir + '/result/final_files'
        if anno_type == 'genic':
            f_load = final_dir + '/sfigs/overlap_{}.txt'.format(annotation)
        else:
            f_load = final_dir + '/sfigs/ccre_overlap_{}.txt'.format(annotation)

        # load total and overlap counts for each annotation
        total_counts = []
        overl_counts = []
        with open(f_load, 'r') as f:
            for line in f:
                tot, ovr = map(int, line.split())
                total_counts.append(tot)
                overl_counts.append(ovr)

        # get fraction of each annotation that is conserved
        fraction_in_cons = []
        for i in xrange(anno_num):
            frc = 1.0 * overl_counts[i] / total_counts[i]
            fraction_in_cons.append(frc)
        frac_in.append(fraction_in_cons)

        # get fraction of conserved falling into each annotation
        n_cons = sum(overl_counts)
        fraction_of_cons = []
        for n in overl_counts:
            fraction_of_cons.append(1.0 * n / n_cons)
        frac_of.append(fraction_of_cons)

    ### PANEL A: plot % of each genic annotation that is conserved
    hgt = 2.16 if anno_type == 'genic' else 2.3
    plt.figure(figsize=(6.5, hgt))
    apos, bpos, cpos = 0.105, 0.435, 0.765
    typos = 0.885
    xtckpad, ytckpad = 0.05, 0.05
    ylblpad = 0
    bot = 0.35 if anno_type == 'genic' else 0.38
    top = 0.95 if anno_type == 'genic' else 0.95
    plt.subplots_adjust(top=top, right=0.995, left=0.1, bottom=bot,
                        wspace=0.4)
    ax1 = plt.subplot(131)
    format_panels(ax1)
    xi = np.arange(anno_num)
    off = -0.4
    for fin in frac_in:
        plt.bar(xi+off, fin, width=0.4)
        off += 0.4

    if anno_type == 'genic':
        xlabs = ['CDS', 'UTR', 'UP &\nDOWN', 'SPLICE', 'INTRON', 'INTERG']
    else:
        xlabs = ['CDS', 'PLS', 'ELS', 'CTCF', 'H3K4me3', 'NONE']

    plt.xticks(xi-0.2, xlabs, y=xtckpad, rotation=90)
    plt.xlabel('annotation')
    plt.yticks(x=ytckpad)
    plt.ylabel('fraction IN 6% most\nconstrained sites', labelpad=ylblpad)
    plt.ylim(0, 0.85)
    plt.legend(['phastCons', 'CADD'])
    plt.text(apos, typos, 'a', transform=plt.gcf().transFigure,
             fontweight='bold')

    ### PANEL B: plot % of conserved in each category
    ax2 = plt.subplot(132)
    format_panels(ax2)
    xi = np.arange(anno_num)
    off = -0.4
    for fof in frac_of:
        plt.bar(xi+off, fof, width=0.4)
        off += 0.4
    if anno_type == 'genic':
        xlabs = ['CDS', 'UTR', 'UP &\nDOWN', 'SPLICE', 'INTRON', 'INTERG']
    else:
        xlabs = ['CDS', 'PLS', 'ELS', 'CTCF', 'H3K4me3', 'NONE']
    plt.xticks(xi-0.2, xlabs, y=xtckpad, rotation=90)
    plt.xlabel('annotation')
    plt.yticks(x=ytckpad)
    plt.ylabel('fraction OF 6% most\nconstrained sites', labelpad=ylblpad)
    plt.ylim(0, 0.75)
    plt.text(bpos, typos, 'b', transform=plt.gcf().transFigure,
             fontweight='bold')

    ### PANEL C: plot enrichment in each category
    ax3 = plt.subplot(133)
    format_panels(ax3)
    xi = np.arange(anno_num)
    off = -0.4
    for fin in frac_in:
        plt.bar(xi + off, [fr/0.06 for fr in fin], width=0.4)
        off += 0.4
    plt.axhline(1, ls='--', lw=0.5, color='k')
    if anno_type == 'genic':
        xlabs = ['CDS', 'UTR', 'UP &\nDOWN', 'SPLICE', 'INTRON', 'INTERG']
    else:
        xlabs = ['CDS', 'PLS', 'ELS', 'CTCF', 'H3K4me3', 'NONE']
    plt.xticks(xi-0.2, xlabs, y=xtckpad, rotation=90)
    plt.xlabel('annotation')
    plt.yticks(x=ytckpad)
    plt.ylabel('fold-enrichment in 6%\nmost constrained sites',
               labelpad=ylblpad)
    plt.legend(['1.0'])

    if anno_type == 'genic':
        f_save = final_dir + '/sfigs/fig_S18.enrichment_of_CADDcons_genic.png'
    else:
        f_save = final_dir + '/sfigs/fig_S18.enrichment_of_CADDcons_cCRE.png'

    plt.text(cpos, typos, 'c', transform=plt.gcf().transFigure,
             fontweight='bold')

    plt.savefig(f_save, dpi=512)
    plt.close()


# plot_combine_cadd_cons()
plot_combine_cadd_cons('ccre')


#%% count numbers of chromhmm annotations per base
cell_lines = ['GM12878', 'H1hesc', 'Huvec', 'Hepg2', 'K562']
states = ['txn', 'pro', 'enh', 'ins', 'rep']
cst = ChromStruct('chr1')
counter_list = []
# for st in states:
state_counter = Counter()
state_arrays = []
for st in states:
    chrom_arrays = []
    for ch in human_autosomes[:1]:
        cst.chrom = ch
        carr = np.zeros(shape=cst.chlen, dtype=int)
        for cl in cell_lines:
            anno = cl.lower() + '_' + st
            for (start, end) in np.loadtxt(cst.bs_target(anno), usecols=(1,2)).astype(int):
                carr[start:end] += 1
        chrom_arrays.append(carr)
    sarr = np.concatenate(chrom_arrays)
    state_arrays.append(sarr)
            # state_counter += Counter(carr)
    # counter_list.append(state_counter)


#%%
bcounts = [np.bincount(ar) for ar in state_arrays]
bcnt = [b[1:] for b in bcounts]
#%%
labels = ['transcription', 'promoter', 'enhancer']
for i in xrange(3):
    b = bcnt[i]
    label = states[i]
    x = np.arange(len(b)) + 1
    plt.figure(figsize=(3.25, 2.5))
    plt.subplots_adjust(top=0.9, right=0.99, left=0.19, bottom=0.2)
    plt.bar(x, b / (1.0*np.sum(b)))
    plt.title(labels[i])
    plt.xlabel('number of overlapping cell lines')
    plt.ylabel('percent of annotation')
    final_dir = root_dir + '/result/final_files'
    f_save = final_dir + '/sfigs/{}.counts.png'.format(label)
    plt.savefig(f_save, dpi=512)
    plt.close()


#%% make an hg38 chrom size dict
hg38_file = root_dir+ '/data/coords/hg38.chrom.sizes'
hg38_size = {}
with open(hg38_file, 'r') as f:
    for line in f:
        chrom, size = line.split()
        if chrom in human_autosomes:
            hg38_size[chrom] = int(size)


#%% make dict of hg38 length empty arrays for each chrom
anno_dict = {}
for anno in ['PLS', 'ELS', 'CTCF']:
    chrom_dict = {}
    for ch in human_autosomes:
        ch_arr = np.zeros(shape=hg38_size[ch], dtype='u2')
        chrom_dict[ch] = ch_arr
    anno_dict[anno] = chrom_dict

# go through each cell-specific cCRE file and fill out each chrom
ccre_folder = root_dir + '/data/coords/cCRE_by_type/'
for fccre in os.listdir(ccre_folder):
    if fccre.endswith('.bed'):
        fccre = ccre_folder + fccre
        with open(fccre, 'r') as f:
            for line in f:
                line = line.split()
                chrom = line[0]
                if chrom not in human_autosomes:
                    continue
                anno = line[-2]
                start, end = map(int, line[1:3])
                primary_anno = anno.split(',')[0]
                if primary_anno in ['dELS', 'pELS']:
                    anno_dict['ELS'][chrom][start:end] += 1
                elif primary_anno == 'PLS':
                    anno_dict['PLS'][chrom][start:end] += 1
                elif primary_anno == 'CTCF-only':
                    anno_dict['CTCF'][chrom][start:end] += 1
                else:
                    continue



#%%
annos = ['PLS', 'ELS', 'CTCF']
overlap_count = {}
for anno in annos:
    cellcount = np.zeros(shape=25)
    for ch in human_autosomes:
        a = anno_dict[anno][ch]
        for i in xrange(25):
            cellcount[i] += np.sum(a == i+1)
    overlap_count[anno] = cellcount

#%%
plt.figure(figsize=(10,10))
xi = np.arange(25)
for i, anno in enumerate(annos, 1):
    plt.subplot(3, 1, i)
    plt.title(anno)
    plt.bar(xi+1, overlap_count[anno])
plt.show()
#%%
# DNase-H3K4me3 4736
# CTCF-only,CTCF-bound 14094
# DNase-H3K4me3,CTCF-bound 1652
# PLS 12485
# pELS,CTCF-bound 1925
# PLS,CTCF-bound 1957
# Low-DNase 785308
# pELS 15961
# DNase-only 48644
# dELS 36772
# dELS,CTCF-bound 3001

#%% (NOT USED IN PAPER)
def plot_overlap_analysis(annotation, label):
    # set path to saved data
    final_dir = root_dir + '/result/final_files'
    f_load = final_dir + '/sfigs/overlap_{}.txt'.format(annotation)

    # load total and overlap counts for each annotation
    total_counts = []
    overl_counts = []
    with open(f_load, 'r') as f:
        for line in f:
            tot, ovr = map(int, line.split())
            total_counts.append(tot)
            overl_counts.append(ovr)

    # get fraction of each annotation that is conserved
    fraction_in_cons = []
    for i in xrange(6):
        frc = 1.0 * overl_counts[i] / total_counts[i]
        fraction_in_cons.append(frc)

    # get fraction of conserved falling into each annotation
    n_cons = sum(overl_counts)
    fraction_of_cons = []
    for n in overl_counts:
        fraction_of_cons.append(1.0*n/n_cons)

    # plot % of each genic annotation that is conserved
    plt.figure(figsize=(3.25, 2.5))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.19, bottom=0.3)
    ax1 = plt.subplot(111)
    format_panels(ax1)
    xi = range(6)
    plt.bar(xi, fraction_in_cons)
    xlabs = ['CDS', 'UTR', 'UP\nDOWN', 'SPLICE', 'INTRON', 'INTERG']
    plt.xticks(xi, xlabs, y=0.04, rotation=90)
    plt.xlabel('annotation')
    plt.yticks(x=0.03)
    plt.ylabel('fraction IN 6% {}'.format(label))
    plt.ylim(0, 0.85)

    f_save = final_dir + '/sfigs/fraction_genic_in_{}.png'.format(label)
    plt.savefig(f_save, dpi=512)
    plt.close()

    # plot % of conserved in each category
    plt.figure(figsize=(3.25, 2.5))
    plt.subplots_adjust(top=0.99, right=0.99, left=0.19, bottom=0.3)
    ax1 = plt.subplot(111)
    format_panels(ax1)
    xi = range(6)
    plt.bar(xi, fraction_of_cons)
    xlabs = ['CDS', 'UTR', 'UP\nDOWN', 'SPLICE', 'INTRON', 'INTERG']
    plt.xticks(xi, xlabs, y=0.04, rotation=90)
    plt.xlabel('annotation')
    plt.yticks(x=0.03)
    plt.ylabel('fraction OF 6% {}'.format(label))
    plt.ylim(0, 0.55)
    f_save = final_dir + '/sfigs/fraction_of_{}_genic.png'.format(label)
    plt.savefig(f_save, dpi=512)
    plt.close()


a1 = 'fish_cons94_gmask'
a2 = 'cadd94_gmask'
plot_overlap_analysis(a1, 'conserved')
plot_overlap_analysis(a2, 'CADD')

