__author__ = 'davidmurphy'


import os
import numpy as np
from collections import defaultdict
from data_processing.data_tools import binary_mask_segments
from classes.runstruct import root_dir, human_autosomes, chromosome_length
ccre_file = root_dir + '/data/coords/V2.hg19-liftOver-cCREs-Final.bed'


class CCRERecord:
    """object representation of a line of cCRE data from ENCODE Project"""
    def __init__(self, line):
        # assign member variables for each field in the line
        fields = line.split()
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.rdhs_accession = fields[3]
        self.ccre_accession = fields[4]
        self.group = fields[5]

    @property
    def primary_mark(self):
        """get the primary mark for multi-label cCREs like pELS,CTCF-bound"""
        return self.group.split(',')[0]

    @property
    def writeline(self):
        """formated output line for saving data as bs annos in bed format"""
        return '{} {} {}\n'.format(self.chrom, self.start, self.end)


class CCREList:
    """a class containing a list of cCREs of a given type"""
    def __init__(self, group):
        # set group type and create an empty list for storing data
        self.group = group
        self.ccrelist = []

    @property
    def totalsites(self):
        """get the total site count in the cCRE list"""
        return sum(c.end - c.start for c in self.ccrelist)

    @property
    def meanlength(self):
        """get mean length of cCREs in list"""
        return self.totalsites / len(self.ccrelist)


#%% MAKE GROUPS OF COMBINED CCRES
def make_basic_groups():
    """make a basic set of the 3 groups of interest: PLS, ELS, CTCF"""
    # create dictionaries of lists for each cCRE-type for each chrom
    els_dict = dict((ch, CCREList('ELS')) for ch in human_autosomes)
    pls_dict = dict((ch, CCREList('PLS')) for ch in human_autosomes)
    ctcf_dict = dict((ch, CCREList('CTCF')) for ch in human_autosomes)
    me3_dict = dict((ch, CCREList('H3K4me3')) for ch in human_autosomes)

    #DNase-H3K4me3
    # put each line in its respective list and dictionary for the 3 classes
    with open(ccre_file, 'r') as f:
        for line in f:
            # create a CCRERecords object
            ccre = CCRERecord(line)
            # skip wierd liftover non-standard chroms
            if ccre.chrom not in human_autosomes:
                continue
            # enhancer like signatures (distal or proximal)
            if ccre.primary_mark in ['dELS', 'pELS']:
                els_dict[ccre.chrom].ccrelist.append(ccre)
            # promoter like signatures
            elif ccre.primary_mark == 'PLS':
                pls_dict[ccre.chrom].ccrelist.append(ccre)
            # CTCF bound (only)
            elif ccre.primary_mark == 'CTCF-only':
                ctcf_dict[ccre.chrom].ccrelist.append(ccre)
            # DNAse-H3K4me3 marks
            elif ccre.primary_mark == 'DNase-H3K4me3':
                me3_dict[ccre.chrom].ccrelist.append(ccre)
            else:
                continue

    return els_dict, pls_dict, ctcf_dict, me3_dict


#%% save original combined cCREs to coords dir
# cre_lbls = ['ELS', 'PLS', 'CTCF', 'H3K4me3']
# cre_dicts = [els, pls, ctcf, me3]
# # cre_lbls = []
# # cre_dicts = [me3]
# for (l, d) in zip(cre_lbls, cre_dicts):
#     anno = 'cCRE_{}'.format(l)
#     d_out = root_dir + '/data/coords/{}'.format(anno)
#     if not os.path.isdir(d_out):
#         os.mkdir(d_out)
#     for ch in human_autosomes:
#         f_out = d_out + '/{}.{}.bed'.format(ch, anno)
#         with open(f_out, 'w') as f:
#             for ccre in d[ch].ccrelist:
#                 f.write(ccre.writeline)
#
#
# #%% check the degree of overlap within cCREs
# arr_dict = {}
# for ch in human_autosomes:
#     # create zero array of chrom size
#     a = np.zeros(shape=chromosome_length(ch), dtype=int)
#     # go through each cCRE dict
#     for d in [els, pls, ctcf]:
#         # go through each cCRE line for the current chrom and cCRE
#         for ccre in d[ch].ccrelist:
#             a[ccre.start:ccre.end] += 1
#
#     # put the array in the dictionary for the current chrom
#     arr_dict[ch] = a
#
#
# # check the overlap between cCREs and genic annotations
# counts = []
# for ch in human_autosomes:
#     a = np.copy(arr_dict[ch])
#     fcds = root_dir + '/data/bsanno/cds/{}.cds.bed'.format(ch)
#     # fperi = root_dir + '/data/bsanno/peri/{}.peri.bed'.format(ch)
#     # for f in [fcds, fperi]:
#     f = fcds
#     for (start, end) in np.loadtxt(f, usecols=(1,2), dtype=int):
#         a[start:end] += 1
#
#     # trim array to nonzeros only
#     b = a[a>0]
#     counts.append(b)
#
# b = np.concatenate(counts)
# h = np.histogram(b, bins=[0,1,2,3,4])
#
# #%%
# els_sum = sum(ccrelist.totalsites for ccrelist in els.values())
# pls_sum = sum(ccrelist.totalsites for ccrelist in pls.values())
# ctcf_sum = sum(ccrelist.totalsites for ccrelist in ctcf.values())


#%% make a set of annotations after removing overlaps with genic/splice
def genic_filtered_ccre(chrom, cre_lbl):
    # set path to UN-filtered cCRE file
    d_coords = root_dir + '/data/coords'
    d_in = d_coords + '/cCRE_{}'.format(cre_lbl)
    f_in = d_in + '/{}.cCRE_{}.bed'.format(chrom, cre_lbl)

    # convert annotation to a chrom length logical array
    a = np.zeros(shape=chromosome_length(chrom), dtype=bool)

    # flip all sites in the UN-filtered data array to 1
    for (start, end) in np.loadtxt(f_in, usecols=(1,2), dtype=int):
        a[start:end] = True

    # calculate the initial unfiltered coverage of chrom
    init_cov = np.sum(a)

    # now load cds and peri annotations to flip sites back to 0s that overlap
    # for an in ['cds', 'peri']:
    for an in ['cds']:
        f_in = root_dir + '/data/bsanno/{a}/{c}.{a}.bed'.format(a=an, c=chrom)
        for (start, end) in np.loadtxt(f_in, usecols=(1,2), dtype=int):
            a[start:end] = False

    # calculate the filtered coverage of the chrom
    end_cov = np.sum(a)

    # convert the mask to segments
    segs = binary_mask_segments(a)

    # set path to filtered cCRE file that will be written
    anno = 'cCRE_{}_filtered'.format(cre_lbl)
    d_out = root_dir + '/data/bsanno/{}'.format(anno)
    if not os.path.isdir(d_out):
        os.mkdir(d_out)
    f_out = d_out + '/{}.{}.bed'.format(chrom, anno)

    # write new annotation to file
    with open(f_out, 'w') as f:
        for (start, end) in segs:
            f.write('{} {} {}\n'.format(chrom, start, end))

    # print status message showing the fraction removed by filters
    msg = '{} {}: initial={} final={} fraction_removed={}'
    frac = 1.0 * end_cov / init_cov
    print msg.format(chrom, cre_lbl, init_cov, end_cov, 1.0 - frac)

    return None


#%%
# cre_lbls = ['ELS', 'PLS', 'CTCF', 'H3K4me3']
# for l in cre_lbls:
#     for ch in human_autosomes:
#         genic_filtered_ccre(ch, l)
#     print '----'
#

#%% SPLIT ANNOTATIONS BY DEGREE OF OVERLAP
def get_cre_label(lbl):
    """function to get labels from cCRE file lines"""
    full_label = lbl.split(',')[0]
    if full_label in ['dELS', 'pELS']:
        return 'ELS'
    # promoter like signatures
    elif full_label == 'PLS':
        return 'PLS'
    # CTCF bound (only)
    elif full_label == 'CTCF-only':
        return 'CTCF'
    # DNAse-H3K4me3 marks
    elif full_label == 'DNase-H3K4me3':
        return 'H3K4me3'
    else:
        print 'ERROR: label {} is not recognized'.format(full_label)


def cds_filter(chrom):
    """flip cds overlapping sites to 0s for a chrom array"""
    a = np.ones(shape=chromosome_length(chrom), dtype=bool)
    f_in = root_dir + '/data/bsanno/{a}/{c}.{a}.bed'.format(a='cds', c=chrom)
    for (start, end) in np.loadtxt(f_in, usecols=(1, 2), dtype=int):
        a[start:end] = False
    return a


def save_segments(chrom, segs, f_out):
    # write new annotation to file
    with open(f_out, 'w') as f:
        for (start, end) in segs:
            f.write('{} {} {}\n'.format(chrom, start, end))


# def split_ccre_by_overlap_count()
# for each anno, make a chromosome array dict
cre_lbls = ['ELS', 'PLS', 'CTCF', 'H3K4me3']
anno_dict = {}
for anno in cre_lbls:
    chrom_dict = {}
    # for each chromosome, create an empty integer array of chrom length
    for ch in human_autosomes:
        ch_arr = np.zeros(shape=chromosome_length(ch), dtype='u1')
        # put the empty chromosome array into the chrom_dict
        chrom_dict[ch] = ch_arr
    # put the dict of empty arrays for the current anno into anno_dict
    anno_dict[anno] = chrom_dict


# create arrays of overlap counts for cCREs
# iterate over cell-type specific cCREs and add them to arrays
ccre_dir = root_dir + '/data/coords/cCRE_by_type'
f_ccre = ccre_dir + '/cCREs.combined.across.25.cell.types.liftOver.to.hg19.bed'
with open(f_ccre, 'r') as f:
    for line in f:
        line = line.split()
        # get the label for cCRE type
        anno = get_cre_label(line[3])
        # get chrom and coords
        ch = line[0]
        # skipped unmapped chroms
        if ch not in human_autosomes:
            continue
        start, end = map(int, line[1:3])
        # use anno, chrom and coords to add info across cell types to array
        anno_dict[anno][ch][start:end] += 1


#%% save annotations to bsanno dir
# for each anno, go through chroms and create annos where counts are above and
# below the medians
median_dict = {'PLS':14, 'ELS':1, 'H3K4me3':1, 'CTCF':3}
bsanno_dir = root_dir + '/data/bsanno'
for anno in cre_lbls:
    for ch in human_autosomes:
        # get chrom specific filter array
        filter_arr = cds_filter(ch)
        # get anno and chrom specific array
        arr = anno_dict[anno][ch]

        # # get array of sites where the counts are less than or equal to median
        # # filter overlapping cds regions
        # below_arr = ((0 < arr) & (arr <= median_dict[anno]) & filter_arr)
        # # convert array to segments
        # below_segs = binary_mask_segments(below_arr)
        # # create file save paths and file names
        # below_dir = bsanno_dir + '/{}_below_median'.format(anno)
        # if not os.path.isdir(below_dir):
        #     os.mkdir(below_dir)
        # below_file = below_dir + '/{}.{}_below_median.bed'.format(ch, anno)
        # save_segments(ch, below_segs, below_file)
        #
        # # get array of sites where the counts are greater than median
        # # filter overlapping cds regions
        # above_arr = ((arr > median_dict[anno]) & filter_arr)
        # above_segs = binary_mask_segments(above_arr)
        # # create file save paths and file names
        # above_dir = bsanno_dir + '/{}_above_median'.format(anno)
        # if not os.path.isdir(above_dir):
        #     os.mkdir(above_dir)
        # above_file = above_dir + '/{}.{}_above_median.bed'.format(ch, anno)
        # save_segments(ch, above_segs, above_file)

        # create filtered integrated annotation
        union_arr = ((0 < arr) & filter_arr)
        union_segs = binary_mask_segments(union_arr)
        union_dir = bsanno_dir + '/{}_union'.format(anno)
        if not os.path.isdir(union_dir):
            os.mkdir(union_dir)
        union_file = union_dir + '/{}.{}_union.bed'.format(ch, anno)
        save_segments(ch, union_segs, union_file)


#%% count annotation types
tot_sites = 0
for anno in cre_lbls:
    n_sites = 0
    for ch in human_autosomes:
        # union_dir = bsanno_dir + '/{}_union'.format(anno)
        # union_file = union_dir + '/{}.{}_union.bed'.format(ch, anno)
        # n_union = 0
        # for start, end in np.loadtxt(union_file, usecols=(1,2)):
        #     n_union += (end-start)
        # n_sites += n_union
        # msg = '{} {} union={}'
        # print msg.format(anno, ch, n_union)
        below_dir = bsanno_dir + '/{}_below_median'.format(anno)
        below_file = below_dir + '/{}.{}_below_median.bed'.format(ch, anno)
        n_below = 0
        for start, end in np.loadtxt(below_file, usecols=(1,2)):
            n_below += (end-start)
        above_dir = bsanno_dir + '/{}_above_median'.format(anno)
        above_file = above_dir + '/{}.{}_above_median.bed'.format(ch, anno)
        n_above = 0
        for start, end in np.loadtxt(above_file, usecols=(1, 2)):
            n_above += (end - start)
        n_sites += (n_below+n_above)

        msg = '{} {} below_median={} above_median={}'
        print msg.format(anno, ch, n_below, n_above)
    msg = '{} total_sites={} fraction_hg19={:.2f}%'
    print msg.format(anno, n_sites, n_sites/2.88e7)
    tot_sites += n_sites
msg = 'all_cCRE total_sites={} fraction_hg19={:.2f}%'
print msg.format(tot_sites, tot_sites/2.88e7)

#%%
overlap_count = {}
for anno in cre_lbls:
    cellcount = np.zeros(shape=25)
    for ch in human_autosomes:
        a = anno_dict[anno][ch]
        for i in xrange(25):
            cellcount[i] += np.sum(a == i+1)
    overlap_count[anno] = cellcount


#%%
import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
plt.subplots_adjust(hspace=0.3)
xi = np.arange(25)
lbls = 'PLS ELS H3K4me3 CTCF'.split()
for i, anno in enumerate(lbls, 1):
    acnt = np.array(overlap_count[anno])
    acnt /= 283840.0
    plt.subplot(4, 1, i)
    plt.title(anno)
    plt.bar(xi+1, acnt)
final_dir = root_dir + '/result/final_files'
f_save = final_dir + '/sfigs/cCRE.counts.png'
plt.savefig(f_save, dpi=512)


#%%
# # count sites in CADD exon/nonexon
# caddex = 30748862.0
# caddnex = 142210548.0
# cadd_tot = caddex+caddnex
# frac_caddex = caddex / cadd_tot
# frac_caddnex = 1-frac_caddex
# print 'fraction of CADD 6% sites: exonic={:.3f} nonexonic={:.3f}'.format(frac_caddex, frac_caddnex)
#
# # count mutations in CADD exon/nonexon
# cadd_u = [6.847319059223717e-09, 8.990076634381303e-09]
# cadd_uex = caddex * cadd_u[0]
# cadd_unex = caddnex * cadd_u[1]
# cadd_utot = cadd_uex+cadd_unex
# print 'fraction of CADD 6% U: exonic={:.3f} nonexonic={:.3f}'.format(cadd_uex/cadd_utot, cadd_unex/cadd_utot)
#
# # count sites in conserved exon/nonexon
# conex = 26831935.0
# connex = 146111596.0
# con_tot = conex+connex
# frac_conex = conex / con_tot
# frac_connex = 1-frac_conex
# print 'fraction of vertebrate conserved 6% sites: exonic={:.3f} nonexonic={:.3f}'.format(frac_conex, frac_connex)
#
# # count mutations in conserved exon/nonexon
# con_u = [1.4385591476919668e-08, 8.27817203389766e-09]
# con_uex = conex * con_u[0]
# con_unex = connex * con_u[1]
# con_utot = con_uex+con_unex
# print 'fraction of vertebrate conserved 6% U: exonic={:.3f} nonexonic={:.3f}'.format(con_uex/con_utot, con_unex/con_utot)

#%%



