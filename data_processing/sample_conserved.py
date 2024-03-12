__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv
from collections import Counter
from classes.runstruct import ChromStruct, RunStruct
from data_processing.data_tools import mask_segments, binary_mask_segments, \
    conservation_mask, dist_percentile, wigz_array
from classes.geneticmap import GeneticMap


def downsample_conserved(ch, idx):
    """downsample conserved data and save to new files with new init file"""
    # create suffix for new file and annotation labels
    suf = '{:02}'.format(idx)

    # set old/new token labels
    old_tkn = 'pr95.cleanrun'
    new_tkn = old_tkn + '.ds.' + suf

    # set old/new
    b_anno = 'primate_cons95_Segments'
    new_b_anno = b_anno + '_ds_' + suf

    # load ChromStruct and new Chromstruct to save later
    cst = ChromStruct(chrom=ch, tkn=old_tkn)
    rst = RunStruct(tkn=new_tkn, bs_annos=(new_b_anno,), bdir=new_b_anno)

    # set old/new anno file names
    anno_file = cst.bs_target(b_anno)
    new_anno_file = anno_file.replace(b_anno, new_b_anno)


    # make new dir, if needed
    new_dir = '/'.join(new_anno_file.split('/')[:-1])
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)

    # load conserved segments
    cons = np.loadtxt(anno_file, dtype=str)

    # create downsample mask
    msize = len(cons)
    pq = [0.9, 0.1]
    choices = [True, False]
    smpl = np.random.choice(choices, size=msize, p=pq)

    # remove masked conserved sites with downsample mask
    mcons = cons[smpl]

    # save downsampled conserved sites to new file
    np.savetxt(new_anno_file, mcons, fmt='%s %s %s')

    # save indices for easy use later
    idx_file = new_anno_file.replace('.bed', '.sidx.npy')
    np.save(idx_file, smpl)

    # save init file if it does not exist yet
    if not os.path.isfile(rst.txt_file):
        rst.save()

    return None


def filter_neutral_conserved(ch, b_anno):
    """remove conserved sites marked as neutral in neutmask"""
    # create suffix for new file and annotation labels
    suf = 'neut'

    # # set old/new token labels
    # old_tkn = 'pr95.nff_2'
    # # new_tkn = old_tkn + '.rm.' + suf
    # new_tkn = 'pr92.nff2.rmneut'

    # set old/new
    new_b_anno = b_anno + '_rm_' + suf

    # load ChromStruct and new Chromstruct to save later
    # cst = ChromStruct(chrom=ch, tkn=old_tkn)
    cst = ChromStruct(chrom=ch, tkn='basic.ncons')

    # set old/new anno file names
    anno_file = cst.bs_target(b_anno)
    new_anno_file = anno_file.replace(b_anno, new_b_anno)

    # make new dir, if needed
    new_dir = '/'.join(new_anno_file.split('/')[:-1])
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)

    # load conserved segments as int
    cons = np.loadtxt(anno_file, usecols=(1, 2), dtype=int)

    # convert conserved segments into a conserved mask with 1s for cons
    cons_mask = mask_segments(np.zeros(shape=cst.chlen), cons, flipoff=False)
    n_1 = cons_mask.sum()

    # get index of neutral positions for chromosome from neutral mask
    nidx = np.where(np.load(cst.neut_masks)['neutmask'] > 0)[0]

    # set neutral positions to 0s and convert back to segments
    cons_mask[nidx] = 0
    n_2 = cons_mask.sum()
    cons_segs = binary_mask_segments(cons_mask)

    # print change in conserved segments
    dpct = 100.0 * (n_1 - n_2) / n_1
    print '{} reduction in conserved sites = {:.3f}%'.format(ch, dpct)

    # create list of new conserved segments
    with open(new_anno_file, 'w') as wf:
        for (start, end) in cons_segs:
            line = '{} {} {}\n'.format(cst.chrom, start, end)
            wf.write(line)

    # # save init file if it does not exist yet
    # rst = RunStruct(tkn=new_tkn, bs_annos=(new_b_anno,), bdir=new_b_anno)
    # if not os.path.isfile(rst.txt_file):
    #     rst.save()

    return None


def find_sitematched_anno(anno):
    """find the annotation that matches cons95 unfiltered most closely"""
    old_tkn = 'pr95.nff_2'
    b_anno = 'primate_cons95_Segments'
    cst = ChromStruct(chrom='chr1', tkn=old_tkn)
    cons95_count = 0
    nmsk = []
    cmsk = []
    for ch in cst.chroms:
        cst.chrom = ch
        # 1. count segment lengths for all 95% cons
        anno_file = cst.bs_target(b_anno)
        for (i, j) in np.loadtxt(anno_file, usecols=(1, 2), dtype=int):
            cons95_count += (j-i)

        # 2. load neutral mask concatenated across autosomes
        nmsk.append(np.load(cst.neut_masks)['neutmask'] > 0)

        # 3. load segments from anno and concatenate as segment mask
        anno_file = cst.bs_target(anno)
        cons = np.loadtxt(anno_file, usecols=(1, 2), dtype=int)
        cms = np.zeros(shape=cst.chlen)
        cmsk.append(mask_segments(cms, cons, flipoff=False))

    # concatenate masks
    nmsk = np.concatenate(nmsk)
    cmsk = np.concatenate(cmsk)

    # check size match and print counts for masks and cons95
    assert cmsk.size == nmsk.size
    print '{} conserved site count: {}'.format(b_anno, cons95_count)
    print '{} neutral site count: {}'.format(old_tkn, nmsk.sum())
    print '{} conserved site count: {}'.format(anno, cmsk.sum())

    # measure the remaining count after filtering neutral sites from anno
    filtered_count = cmsk[~nmsk].sum()
    fraction_95 = 1.0 * filtered_count / cons95_count
    msg = 'fraction of {} in filtered {} segments = {}'
    print msg.format(b_anno, anno, fraction_95)

    return None


def build_conserved_segments(ch, cons, pct):
    """create BED file of segments from the top X-percent conserved sites"""
    # initialize a ChromStruct to retrieve file names for function calls
    cst = ChromStruct(chrom=ch, cons=cons)

    # set paths to conservation scores
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)
    fdist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=cons)

    # get minimum phastCons value for the given percentile
    pmin = dist_percentile(fdist, pct)

    # get conserved status mask for values in the selected percentile
    cmsk = conservation_mask(fwigz, ch, pmin, 1.0, 'u1')

    # turn the conserved site mask into conserved segments
    segs = binary_mask_segments(cmsk)

    # create output path (if needed)
    anno = '{}_cons{}_Segments'.format(cons, int(100*pct))
    fpth = '{}/bsanno/{}'.format(cst.data, anno)
    if not os.path.isdir(fpth):
        os.mkdir(fpth)

    # save new segments
    fout = '{pth}/{c}.{an}.bed'.format(pth=fpth, c=ch, an=anno)
    fmt = '{c} {start} {end}\n'
    with open(fout, 'w') as f:
        for (i, j) in segs:
            f.write(fmt.format(c=ch, start=i, end=j))


def get_filtered_cmask(ch, cons, pct, ncons, npct):
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    # path to conserved distribution and wigz file
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)
    # fdist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=cons)
    fdist = '{p}/{cn}/{cn}.filtered.scores.dist.txt'.format(p=pth, cn=cons)
    # path to conserved distribution and neutral wigz file
    fnwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=ncons, c=ch)
    fndist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=ncons)

    # get maximum neutral score value for neutral percentile
    pmax = dist_percentile(fndist, npct)
    # get neutral mask at the cutoff
    nmsk = conservation_mask(fnwigz, ch, 0.0, pmax)

    # get minimum conserved value for the given percentile
    pmin = dist_percentile(fdist, pct)
    # get conserved status mask for values in the selected percentile
    cmsk = conservation_mask(fwigz, ch, pmin, 1.0)

    # take the intersection of conserved and NOT neutral
    filtered_cmsk = (cmsk & ~nmsk).astype('u1')

    return filtered_cmsk


def build_filtered_conserved_segments(ch, cons, pct, ncons, npct):
    """create BED file of segments from the top X-percent conserved sites"""
    # initialize a ChromStruct to retrieve file names for function calls
    cst = ChromStruct(chrom=ch, cons=cons)

    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    # path to conserved distribution and wigz file
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)
    # fdist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=cons)
    fdist = '{p}/{cn}/{cn}.filtered.scores.dist.txt'.format(p=pth, cn=cons)
    # path to conserved distribution and neutral wigz file
    fnwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=ncons, c=ch)
    fndist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=ncons)

    # get maximum neutral score value for neutral percentile
    pmax = dist_percentile(fndist, npct)
    # get neutral mask at the cutoff
    nmsk = conservation_mask(fnwigz, ch, 0.0, pmax)

    # get minimum conserved value for the given percentile
    pmin = dist_percentile(fdist, pct)
    # get conserved status mask for values in the selected percentile
    cmsk = conservation_mask(fwigz, ch, pmin, 1.0)

    # take the intersection of conserved and NOT neutral
    filtered_cmsk = (cmsk & ~nmsk).astype('u1')

    # print difference from filtering
    print 'initial/filtered: {} {}'.format(cmsk.sum(), filtered_cmsk.sum())

    # turn the conserved site mask into conserved segments
    segs = binary_mask_segments(filtered_cmsk)

    # create output path (if needed)
    fmt = '{}_cons{}_{}_neut{}_filtered'
    anno = fmt.format(cons, int(100*pct), ncons, int(100*npct))
    fpth = '{}/bsanno/{}'.format(cst.data, anno)
    if not os.path.isdir(fpth):
        os.mkdir(fpth)

    # save new segments
    fout = '{pth}/{c}.{an}.bed'.format(pth=fpth, c=ch, an=anno)
    fmt = '{c} {start} {end}\n'
    with open(fout, 'w') as f:
        for (i, j) in segs:
            f.write(fmt.format(c=ch, start=i, end=j))


def build_filtered_conserved_range_segments(ch, cons, pmin, pmax, ncons, npct):
    """create BED file of segments from the top X-percent conserved sites"""
    # initialize a ChromStruct to retrieve file names for function calls
    cst = ChromStruct(chrom=ch, cons=cons)

    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    # path to conserved distribution and wigz file
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)
    # fdist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=cons)
    fdist = '{p}/{cn}/{cn}.filtered.scores.dist.txt'.format(p=pth, cn=cons)
    # path to conserved distribution and neutral wigz file
    fnwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=ncons, c=ch)
    fndist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=ncons)

    # get maximum neutral score value for neutral percentile
    npmax = dist_percentile(fndist, npct)
    # get neutral mask at the cutoff
    nmsk = conservation_mask(fnwigz, ch, 0.0, npmax)

    # get minimum conserved value for the given percentile
    cpmin = dist_percentile(fdist, pmin)
    cpmax = dist_percentile(fdist, pmax)
    # get conserved status mask for values in the selected percentile
    cmsk = conservation_mask(fwigz, ch, cpmin, cpmax)

    # take the intersection of conserved and NOT neutral
    filtered_cmsk = (cmsk & ~nmsk).astype('u1')

    # print difference from filtering
    print 'initial/filtered: {} {}'.format(cmsk.sum(), filtered_cmsk.sum())

    # turn the conserved site mask into conserved segments
    segs = binary_mask_segments(filtered_cmsk)

    # create output path (if needed)
    fmt = '{}_cons{}to{}_{}_neut{}_filtered'
    anno = fmt.format(cons, int(100*pmin), int(100*pmax), ncons, int(100*npct))
    fpth = '{}/bsanno/{}'.format(cst.data, anno)
    if not os.path.isdir(fpth):
        os.mkdir(fpth)

    # save new segments
    fout = '{pth}/{c}.{an}.bed'.format(pth=fpth, c=ch, an=anno)
    fmt = '{c} {start} {end}\n'
    with open(fout, 'w') as f:
        for (i, j) in segs:
            f.write(fmt.format(c=ch, start=i, end=j))


def fish_zeros_dist(ch, cons):
    """count scores at sites where fish PHASTCON is 0"""
    # set file paths
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_fish = '{p}/fish/fish.{c}.wig.gz'.format(p=pth, c=ch)
    f_cons = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)

    # get neutral mask for fish 0s
    fish_mask = conservation_mask(f_fish, ch, 0.0, 0.0)
    # get the scores for "cons" across chromosome
    cons_mask = wigz_array(ch, f_cons)

    # count up each score in the fish 0s
    score_dict = Counter(cons_mask[fish_mask])

    # save scores for the current chromosome
    fmt = '{p}/{c}.{cn}.fish.zeros.dist.txt'
    fout = fmt.format(p=pth, c=ch, cn=cons)
    with open(fout, 'w') as f:
        for k in sorted(score_dict.keys()):
            line = '{} {}\n'.format(k, score_dict[k])
            f.write(line)


def conserved_segment_stats(cons, pct):
    """measure the number of conserved sites before and after filters"""
    # init chromstruct for file paths
    cst = ChromStruct('chr1', cons=cons, )

    # set annotatios for original and neutral filtered files
    b_anno_1 = '{}_cons{}_Segments'.format(cons, pct)
    b_anno_2 = b_anno_1 + '_rm_neut'

    # create counters for segment lengths for each file
    cnt_1, cnt_2 = Counter(), Counter()

    # load each annotation file and record segment lengths
    for ch in cst.chroms:
        cst.chrom = ch
        i1, j1 = np.loadtxt(cst.bs_target(b_anno_1), usecols=(1, 2)).T
        i2, j2 = np.loadtxt(cst.bs_target(b_anno_2), usecols=(1, 2)).T
        # use counters to summarize segment lengths
        cnt_1 += Counter(j1-i1)
        cnt_2 += Counter(j2-i2)

    # write counts for segment lengths to file
    pth = '/'.join(cst.bs_target(b_anno_1).split('/')[:-2])
    fout = '{}/{}.stats.txt'.format(pth, b_anno_1)

    all_keys = cnt_1.keys() + cnt_2.keys()
    all_keys.sort()
    last_key = -1
    with open(fout, 'w') as f:
        f.write('#size original filtered\n')
        for k in all_keys:
            if k != last_key:
                last_key = k
                line = '{} {} {}\n'.format(int(k) , cnt_1[k], cnt_2[k])
                f.write(line)


def filtered_conserved_distr(ch, cons, ncons, npct):
    """get the distribution of conservation scores after filtering neut"""
    # exit with error if cons == ncons
    if cons == ncons:
        print 'ERROR: CANNOT FILTER {} FROM {}'.format(cons, ncons)
        exit(1)

    # set file paths
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    # path to target phastCons scores to be filtered
    fwigz = '{p}/{cn}/{cn}.{c}.wig.gz'.format(p=pth, cn=cons, c=ch)
    # path phastCons scores used to standardize target score disctribution
    fnwigz = '{p}/{ncn}/{ncn}.{c}.wig.gz'.format(p=pth, ncn=ncons, c=ch)
    # path to score distribution for standardized scores
    fndist = '{p}/{cn}/{cn}.scores.dist.txt'.format(p=pth, cn=ncons)

    # get maximum neutral score value for neutral percentile
    pmax = dist_percentile(fndist, npct)
    # get neutral mask at the cutoff
    nmsk = conservation_mask(fnwigz, ch, 0.0, pmax)

    # load conservation scores for conserved annotation
    cscores = wigz_array(ch, fwigz)
    # set neutral sites to 0 in conservation scores
    cscores[nmsk > 0] = 0

    # count up each type of score
    score_dict = Counter(cscores)

    # save scores for the current chromosome
    fmt = '{}/{}/{}.{}.{}.filtered.scores.dist.txt'
    fout = fmt.format(pth, cons, ch, ncons, int(100*npct))
    with open(fout, 'w') as f:
        for k in sorted(score_dict.keys()):
            line = '{} {}\n'.format(k, score_dict[k])
            f.write(line)


def summarize_neutral(cons):
    cst_1 = ChromStruct('chr1', tkn='basic.ncons')
    cst_2 = ChromStruct('chr1', tkn='fish95_{}50'.format(cons))
    n1 = 0
    n2 = 0
    for ch in cst_1.chroms:
        cst_1.chrom = ch
        cst_2.chrom = ch
        n1 += np.load(cst_1.neut_masks)['neutmask'].sum()
        n2 += np.load(cst_2.neut_masks)['neutmask'].sum()

    msg = '100-VERT NEUTRAL SUM: {:.4e} {} NEUTRAL SUM: {:.4e}'
    print msg.format(n1, cons, n2)


def check_seg_lengths(fname):
    """sum the segment lengths of conserved segments in a folder and print"""
    # get segment files from folder name
    bdir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/bsanno/'
    fldr = bdir + fname + '/'
    files = [fldr + f for f in os.listdir(fldr) if f.endswith('.bed')]
    # check that there is one file per autosome
    assert len(files) == 22
    # sum the segment lengths
    n = 0
    for f in files:
        start, end = np.loadtxt(f, usecols=(1, 2), dtype=int).T
        n += np.sum(end-start)

    print n


def genmap_bin_cons(ch, pct, bin_size):
    """count conserved sites in genetic map unit bins on different scales"""
    # init ChromStruct and GeneticMap classes
    cst = ChromStruct(ch)
    gmp = GeneticMap(ch, cst.gmap_files)

    # create bins for counting up conserved sites in each segment
    gmin, gmax = gmp.gmap_limits
    bins = np.arange(gmin, gmax + bin_size, bin_size)

    # get conserved segments from each species and store genetic positions
    cnt = []
    species = ['ape', 'euarchontoglires', 'fish']
    for spc in species:
        # load segments for annotation
        anno = '{}_cons{}_euarchontoglires_neut30_filtered'.format(spc, pct)
        f_cons = cst.bs_target(anno)
        # convert segments into array of genetic positions
        pos = []
        for (start, end) in np.loadtxt(f_cons, usecols=(1, 2)):
            pos.append(np.arange(start, end, 1))
        pos = np.concatenate(pos)
        gpos = gmp.interp_gpos(pos)
        # sort segment positions into bins
        cnt.append(np.histogram(gpos, bins=bins)[0])

    # save the counts per bin for each species
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    f_save = '{}/{}.{}.gmap.dist.bin_{:.2e}.txt'.format(pth, ch, pct, bin_size)
    l_fmt = '{} {:.0f} {:.0f} {:.0f}\n'
    with open(f_save, 'w') as f:
        for i in range(len(bins)-1):
            l_a, l_e, l_f = cnt[0][i], cnt[1][i], cnt[2][i]
            if all(a == 0 for a in (l_a, l_e, l_f)):
                continue
            else:
                line = l_fmt.format(bins[i], l_a, l_e, l_f)
                f.write(line)


def local_main():
    ch = 'chr21'
    cons = 'primate'
    # filtered_conserved_distr(ch, cons)


def main():
    # if len(argv) != 3:
    #     print 'usage: sample_conserved <chrom> <cons>'
    #     exit(1)
    #
    # ch, cons = argv[1:3]
    # fish_zeros_dist(ch, cons)

    # if len(argv) != 4:
    #     print 'usage: sample_conserved <chrom> <pct> <bin_size>'
    #     exit(1)
    #
    # ch, pct = argv[1:3]
    # bin_size = float(argv[3])
    # genmap_bin_cons(ch, pct, bin_size)

    # if len(argv) != 4:
    #     print 'usage: sample_conserved <chrom> <cons> <pct>'
    #     exit(1)

    # ch, cons = argv[1:3]
    # pct = float(argv[3])
    # ncons = 'euarchontoglires'
    # npct = 0.35
    # # build new filtered and unfiltered segments
    # # build_conserved_segments(ch, cons, pct)
    # if cons != ncons:
    #     build_filtered_conserved_segments(ch, cons, pct, ncons, npct)
    #
    if len(argv) != 5:
        print 'usage: sample_conserved <chrom> <cons> <pmin> <pmax>'
        exit(1)

    ch, cons = argv[1:3]
    pmin = float(argv[3])
    pmax = float(argv[4])
    ncons = 'euarchontoglires'
    npct = 0.35
    # build ranged conserved segs
    build_filtered_conserved_range_segments(ch, cons, pmin, pmax, ncons, npct)

    # if len(argv) != 3:
    #     print 'usage: sample_conserved <chrom> <cons>'
    #     exit(1)
    # ch, cons = argv[1:3]
    # ncons = 'euarchontoglires'
    # npct = 0.30
    # filtered_conserved_distr(ch, cons, ncons, npct)

    # if len(argv) != 4:
    #     print 'usage: sample_conserved <chrom> <cons> <pct>'
    #     exit(1)
    # ch, cons = argv[1:3]
    # pct = float(argv[3])
    # build_conserved_segments(ch, cons, pct)
    # conserved_segment_stats(cons, int(100*pct))
    # filtered_conserved_distr(ch, cons)
    # summarize_neutral(cons)
    #
    # if len(argv) != 2:
    #     print 'usage: sample_conserved <folder_name>'
    #     exit(1)
    # folder = argv[1]
    # check_seg_lengths(folder)


if __name__ == '__main__':
    if os.getcwd().startswith('/ifs/'):
        main()
    else:
        local_main()
