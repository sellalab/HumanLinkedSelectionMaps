__author__ = 'davidmurphy'


import os
import tarfile
import numpy as np
from sys import argv, stderr, stdout
from classes.geneticmap import GeneticMap
from data_processing.precalc_tools import compress_data
from classes.knowngene import complement_strand, complement_base
from classes.runstruct import ChromStruct, root_dir, cst_from_fldr, FileStruct
from data_processing.data_tools import fasta_array, wigz_array, snpcount, \
    cpg_mask


def load_song_maps(ch, pop_id, indv_id='filtered'):
    """load Yun Song's archaic introgression maps"""
    # set directories
    fdir = root_dir + '/data/snps/Song_introgression_maps'
    f_name = fdir + '/{pid}/{pid}_lax_{ch}.tar.gz'.format(pid=pop_id, ch=ch)

    # open the zipped tarball
    tar = tarfile.open(f_name, 'r:gz')
    pos = []
    prob_list = []
    for member in tar.getmembers():
        # load positions
        if '.pos' in str(member):
            f = tar.extractfile(member)
            for line in f:
                pos.append(int(line))
        # load probabilities from each haplotype
        if indv_id in str(member):
            prob = []
            f = tar.extractfile(member)
            for line in f:
                prob.append(float(line))
            prob_list.append(prob)

    # if specific pop ID was used, should only be two haplotypes
    if indv_id != 'filtered':
        assert len(prob_list) == 1

    # convert lists to arrays, take average probs
    prob_list = np.average(np.array(prob_list), axis=0)
    pos = np.array(pos)

    return pos, prob_list


def get_segments(cst):
    """get the segments associated with the Chromstruct"""
    assert isinstance(cst, ChromStruct)

    # create a set of coordinates for the start and end of segments
    segs = np.load(cst.sg_files)['sg']
    end = np.cumsum(segs)
    start = np.concatenate(([0], end[:-1]))

    return np.column_stack((start, end)).astype(int)


def get_nmsk(cst):
    """get the neutral mask associated with the ChromStruct"""
    assert isinstance(cst, ChromStruct)
    print 'loading neutral mask file {}'.format(cst.neut_masks)
    return np.load(cst.neut_masks)['neutmask']


def get_gc_content(cst, nmsk, segs):
    """get the GC content of the neutral sites in nmsk"""
    assert isinstance(cst, ChromStruct)

    # load the reference chromosome
    # ref = fasta_array(cst.chrom, cst.refg_files)
    ref = fasta_array(cst.chrom, cst.ancs_files)

    # get the GC content at neutral sites for each segment
    gc = []
    for (start, end) in segs:
        cur_msk = nmsk[start:end]
        if not np.sum(cur_msk > 0):
            gc.append(0)
        else:
            cur_ref = ref[start:end]
            cur_neut = cur_ref[cur_msk > 0]
            gc_count = np.sum(np.in1d(cur_neut, ['C', 'G']))
            gc_fract = 1.0 * gc_count / len(cur_neut)
            gc.append(gc_fract)

    return np.array(gc)


def get_gc_in_islands(cst, nmsk, segs):
    """get fraction of neutral sites in CpG islands for each segment"""
    assert isinstance(cst, ChromStruct)
    # load the reference chromosome
    # ref = fasta_array(cst.chrom, cst.refg_files)
    ref = fasta_array(cst.chrom, cst.ancs_files)

    # set empty array to fill with CpG sites
    cpg_isl = np.zeros(shape=cst.chlen, dtype=bool)

    # go through each CpG island and flip array to one on CpG island segments
    f_in = cst.data + '/coords/cpg_islands/{}.cpg.islands.bed'.format(cst.chrom)
    with open(f_in, 'r') as f:
        for line in f:
            start, end = map(int, line.split()[1:3])
            cpg_isl[start:end] = True

    # get the fraction of GC out of all neutral sites in CpG islands
    gc_isl_count = []
    for (start, end) in segs:
        # get the segment of the neutral mask
        cur_msk = nmsk[start:end]
        # if there are no neutral sites in the segment, record a 0
        if not np.sum(cur_msk > 0):
            gc_isl_count.append(0)
            continue

        # get the segment of the reference genome
        cur_ref = ref[start:end]
        # get the neutral subset of the segment
        cur_neut = cur_ref[cur_msk > 0]
        # get a mask of GC sites for the neutral sites
        cur_gc = np.in1d(cur_neut, ['C', 'G'])

        # get the segment of the CpG islands array
        isl_seg = cpg_isl[start:end]
        # get the neutral subset of the CpG islands array
        isl_neut = isl_seg[cur_msk > 0]
        # if there are no neutral island sites in the segment, record a 0
        if not np.sum(isl_neut > 0):
            gc_isl_count.append(0)
            continue

        # check that these masks are all the same size
        assert len(cur_neut) == len(cur_gc) == len(isl_neut)
        # get the count of GC neutral sites in islands
        gc_count = (isl_neut & cur_gc)

        # get the fraction of neutral GC sites in islands out of all in islands
        gc_isl_count.append(1.0 * gc_count.sum() / isl_neut.sum())

    return np.array(gc_isl_count)


def get_cpgs_in_islands(cst, nmsk, segs):
    """get fraction of neutral sites in CpG islands for each segment"""
    assert isinstance(cst, ChromStruct)

    # load all CpG sites as a mask
    cpg = np.zeros(shape=cst.chlen, dtype=bool)
    cpg[cpg_mask(cst.chrom, cst.ancs_files)] = True

    # set empty array to fill with CpG sites
    cpg_isl = np.zeros(shape=cst.chlen, dtype=bool)

    # go through each CpG island and flip array to one on CpG island segments
    f_in = cst.data + '/coords/cpg_islands/{}.cpg.islands.bed'.format(cst.chrom)
    with open(f_in, 'r') as f:
        for line in f:
            start, end = map(int, line.split()[1:3])
            cpg_isl[start:end] = True

    # get the fraction of GC out of all neutral sites in CpG islands
    gc_isl_count = []
    for (start, end) in segs:
        # get the segment of the neutral mask
        cur_msk = nmsk[start:end]
        # if there are no neutral sites in the segment, record a 0
        if not np.sum(cur_msk > 0):
            gc_isl_count.append(0)
            continue
        # get current CpG sites
        cur_cpg = cpg[start:end]
        # get the subset of the CpG sites that are neutral
        cpg_neut = cur_cpg[cur_msk > 0]
        # if there are no CpG neutral sites, record 0
        if not np.sum(cpg_neut > 0):
            gc_isl_count.append(0)
            continue
        # get the segment of the CpG islands array
        isl_seg = cpg_isl[start:end]
        # get the neutral subset of the CpG islands array
        isl_neut = isl_seg[cur_msk > 0]
        # if there are no neutral island sites in the segment, record a 0
        if not np.sum(isl_neut > 0):
            gc_isl_count.append(0)
            continue
        # check that these masks are all the same size
        assert len(cpg_neut) == len(isl_neut)

        # get the count of CpG neutral sites in islands
        cpg_in_island = (isl_neut & cpg_neut)

        # get the fraction of neutral GC sites in islands out of all in islands
        gc_isl_count.append(1.0 * cpg_in_island.sum() / cpg_neut.sum())

    return np.array(gc_isl_count)


def get_chrompos(cst, segs):
    """create an array of chrom and segment bounds for each segment"""
    assert isinstance(cst, ChromStruct)
    ch = float(cst.chrom[3:])
    ch_arr = np.full(len(segs), ch)
    chrom_pos_array = np.column_stack((ch_arr, segs)).astype('u4')

    return chrom_pos_array


def get_phast_scores(cst, nmsk, segs):
    """get the phastcons scores for the neutral sites in nmsk"""
    assert isinstance(cst, ChromStruct)
    pth = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/wigz'
    fstr = '{p}/{cn}/{cn}.{c}.wig.gz'
    fwigz = fstr.format(p=pth, cn=cst.ncon, c=cst.chrom)
    con = wigz_array(cst.chrom, fwigz, use_neg1=True)

    # get the mean conservation score at neutral sites for each segment
    mean_con = []
    for (start, end) in segs:
        cur_msk = nmsk[start:end]
        if not np.any(cur_msk > 0):
            mean_con.append(0)
        else:
            cur_con = con[start:end]
            # m = np.mean(cur_con[cur_con >= thresh])
            m = np.mean(cur_con[cur_con >= 0])  # discard unannotated -1 sites
            mean_con.append(m)

    return np.array(mean_con)


def get_cmmb(cst, segs):
    """get the recombination rate from start to end of each segment"""
    assert isinstance(cst, ChromStruct)

    # initialize a GeneticMap and get genetic distances at segment boundaries
    gmp = GeneticMap(cst.chrom, cst.gmap_files)
    gsegs = gmp.interp_gpos(segs)
    gdists = (gsegs[:,1] - gsegs[:,0]) / (segs[:,1] - segs[:,0])

    return gdists


def get_poly_by_type(cst, nmsk, sarr, use_anc=True):
    """organize SNPs by mutation type"""
    # create save directory for SNP types
    snp_dir = root_dir + '/compress/snp_type/'
    if not os.path.isdir(snp_dir):
        os.mkdir(snp_dir)
    # create format for saved files
    save_fmt = snp_dir + cst.chrom + '.{}.{}.npz'

    # get SNP data included ref/alt bases for chrom
    snppos, rbase, rcount, abase, acount = snpcount(cst.snp_files, True)
    # create the array used for compression
    snps = np.column_stack((snppos, rcount, acount))
    # get anc sequence
    anc = fasta_array(cst.chrom, cst.ancs_files)

    # use ancestor sequence to polarize data
    if use_anc:
        # use a different save file name
        save_fmt = snp_dir + cst.chrom + '.{}.{}.fix_anc.npz'
        # get ancestor allele at each snp pos (adjusted for 0-based index)
        snpanc = anc[snppos-1]
        # mask for sites where anc matches alt
        anc_alt = (snpanc == abase)
        # store the alt where it matches anc and ref where it doesnt
        switched_ref = abase[anc_alt]
        switched_alt = rbase[anc_alt]
        # switch alt>ref for sites where alt matches anc
        rbase[anc_alt] = switched_ref
        abase[anc_alt] = switched_alt
        # write err message for number of swithces done
        msg = 'ancestor correction to reference affected {} out of {} sites.\n'
        stderr.write(msg.format(np.sum(anc_alt), len(anc_alt)))
        stdout.flush()

    # combine SNP pairs into 2-char string array
    pairs = np.core.defchararray.add(rbase, abase)

    # for each SNP type, get compressed poly data
    snp_types = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT']
    for snptype in snp_types:
        # make a copy of the neutral mask for filtering originating base types
        cur_nmsk = np.copy(nmsk)
        # # set the base composition background
        # base = snptype[0]
        # # remove neutral sites that don't match the background
        # base_mask = (anc == base) | (anc == complement_base(base))
        # cur_nmsk[~base_mask] = False
        # get complement SNP type
        complement = complement_strand(snptype)
        # create label for current SNP type
        lbl = '{}_{}'.format(snptype, complement)
        # get mask for the snptype and its complement
        smsk = (pairs == snptype) | (pairs == complement)
        # get polymorphism & divergence summaries compressed into segments
        hom, het = compress_data(cur_nmsk, snps[smsk], sarr)[0][:2]
        # save poly data to array
        f_save = save_fmt.format(cst.tkn, lbl)
        np.savez_compressed(f_save, np.column_stack((hom, het)))
        # calculate pi for current SNP type
        pi = 1.0 * np.sum(het) / np.sum(het+hom)
        # write err msg with pi and snp type
        msg = '{} {} pi = {}\n'.format(cst.chrom, lbl, pi)
        stderr.write(msg)
        stdout.flush()

    return None


def get_nonafr_archaic(cst, nmsk, segs, pop_id):
    """get aggregate archaic introgression scores for non-Africans"""
    assert isinstance(cst, ChromStruct)

    # set ID for song maps. "filtered" == use all
    # pop_id = 'NA06986.' + cst.chrom[3:] + '.0'
    indv_id = 'filtered'
    # pop_id = 'CHBS'

    # get segment midpoints and mean probabilities across individuals
    pos, prob = load_song_maps(cst.chrom, pop_id, indv_id)
    # adjust midpoints to start points
    pos = (pos-250).astype(int)
    assert np.all(pos[1:]-pos[:-1] == 500)

    # fill probs for each 500 bp segment
    probs = np.zeros(shape=cst.chlen)
    for i in range(len(pos)):
        pi = pos[i]
        probs[pi:pi+500] = prob[i]

    # for each segment get the mean probability on neutral sites
    seg_probs = []
    for (start, end) in segs:
        mask_seg = nmsk[start:end]
        prob_seg = probs[start:end]
        # seg_probs.append(np.mean(prob_seg[mask_seg]))
        seg_probs.append(np.sum(prob_seg[mask_seg]))

    return np.array(seg_probs)


def get_afr_archaic(cst, nmsk, segs, pop_id):
    """get aggregate archaic introgression scores for Africans"""
    assert isinstance(cst, ChromStruct)

    # set dirs
    snp_dir = root_dir + '/data/snps'
    archie_dir = snp_dir + '/ArchIE-calls-YRI-MSL'

    # use Sriram file for given pop
    use_file = archie_dir + '/{}-freq-map.bed'.format(pop_id)

    # # files from Sriram
    # yri_file = archie_dir + '/YRI-freq-map.bed'
    # msl_file = archie_dir + '/MSL-freq-map.bed'
    # # files from the Durbin paper
    # mean_file = snp_dir + '/YRI.mean.probs.txt'
    # indv_file = snp_dir + '/YRI.NA18933.mean.probs.txt'
    # lezgin_file = snp_dir + '/S_Lezgin-1.archaic.segments.txt'
    # all_file = snp_dir + '/archaic.segments.txt'
    #
    # use_file = yri_file

    # create empty array for probabilities for each site in chromosome
    probs = np.zeros(shape=cst.chlen)

    # get probabilities for the current chromosome
    with open(use_file, 'r') as f:
        for line in f:
            line = line.split()
            # ch, i, j, p = 'chr'+line[1], int(line[2]), int(line[3]), float(line[-6])
            # if ch == cst.chrom and j>i and line[-1] == 'HMM' and line[0] == 'S_Pathan-1':
            ch, i, j, p = line[0], int(line[1]), int(line[2]), float(line[3])
            if ch == cst.chrom:
                probs[i:j] = p

    # for each segment get the mean probability on neutral sites
    seg_probs = []
    for (start, end) in segs:
        mask_seg = nmsk[start:end]
        prob_seg = probs[start:end]
        seg_probs.append(np.mean(prob_seg[mask_seg]))
        # seg_probs.append(np.sum(prob_seg))

    return np.array(seg_probs)


def main():
    if len(argv) != 3:
        print 'usage: segmented_data <folder> <chrom>'
        exit(1)

    # load pre-saved data
    folder, chrom = argv[1:]
    # cst = ChromStruct(chrom=chrom, init=folder)
    cst = cst_from_fldr(folder, ch=chrom)
    cst.files = FileStruct(cst.dict)
    # cst.reset()

    segs = get_segments(cst)
    sarr = np.load(cst.sg_files)['sg']
    nmsk = get_nmsk(cst)

    # save CpG island data
    isl_dir = root_dir + '/compress/il_arr/'
    if not os.path.isdir(isl_dir):
        os.mkdir(isl_dir)
    # np.savez_compressed(cst.il_files, il=get_gc_in_islands(cst, nmsk, segs))
    np.savez_compressed(cst.il_files, il=get_cpgs_in_islands(cst, nmsk, segs))

    # save SNP type data
    get_poly_by_type(cst, nmsk, sarr, use_anc=True)

    # save positional segment data
    sg_dir = root_dir + '/compress/sg_arr/'
    if not os.path.isdir(sg_dir):
        os.mkdir(sg_dir)
    f_save = sg_dir + '{}.sg_arr.{}.npz'.format(chrom, cst.tkn)
    np.savez_compressed(f_save, sg=get_chrompos(cst, segs))

    # # save archaic introgression data
    # ai_dir = root_dir + '/compress/ai_arr/'
    # if not os.path.isdir(ai_dir):
    #     os.mkdir(ai_dir)
    # for pop_id in ['CHBS', 'CEU']:
    #     f_save = cst.ai_files.replace('ai.', 'ai.{}.'.format(pop_id))
    #     ai_arr = get_nonafr_archaic(cst, nmsk, segs, pop_id)
    #     np.savez_compressed(f_save, ai=ai_arr)
    # for pop_id in ['MSL', 'YRI']:
    #     f_save = cst.ai_files.replace('ai.', 'ai.{}.'.format(pop_id))
    #     ai_arr = get_afr_archaic(cst, nmsk, segs, pop_id)
    #     np.savez_compressed(f_save, ai=ai_arr)
    #
    # save GC segment data
    gc_dir = root_dir + '/compress/gc_arr/'
    if not os.path.isdir(gc_dir):
        os.mkdir(gc_dir)
    np.savez_compressed(cst.gc_files, gc=get_gc_content(cst, nmsk, segs))

    # save conservation segment data
    cn_dir = root_dir + '/compress/cn_arr/'
    if not os.path.isdir(cn_dir):
        os.mkdir(cn_dir)
    np.savez_compressed(cst.cn_files, cn=get_phast_scores(cst, nmsk, segs))

    # save cM/Mb segment data
    cm_dir = root_dir + '/compress/cm_arr/'
    if not os.path.isdir(cm_dir):
        os.mkdir(cm_dir)
    np.savez_compressed(cst.cm_files, cm=get_cmmb(cst, segs))


if __name__ == '__main__':
    main()
