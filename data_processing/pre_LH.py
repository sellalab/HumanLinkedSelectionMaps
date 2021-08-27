from precalc_tools import *
from data_tools import *
from sys import argv


__author__ = 'davidmurphy'


def main_local():
    chrom = 'chr22'
    cst = ChromStruct(chrom=chrom,
                      root=root_dir,
                      gmap='AA_Map_10kb',
                      refg='hg19',
                      mask='strictMask',
                      ncon='100vert',
                      cons='primate',
                      outg='rheMac3',
                      neut='YRI',
                      phase='phase3',
                      tkn='pr95.fixB',
                      focal='nonsyn',
                      rm_annos=default_filter_file,
                      bkgd_scale=100,
                      bs_annos=('primate_cons95_Segments',),
                      bdir='bigB100',
                      cs_annos=(),
                      cmap_dir='cmaps',
                      bs_dfe=(default_fitness_effects,),
                      cs_dfe=(),
                      chroms=human_autosomes,
                      methods=('COBYLA', 'Powell'))

    # generate mask files
    mcons = conservation_mask(cst.wigz_ncons, chrom)
    msubs = substitution_mask(cst.axt_outgroup, chrom)
    mcall = sequence_mask(cst.call_mask, chrom)

    # create and save neutral mask
    nmask = neutrality_mask(mcons, mcall, msubs, filter_files=cst.filter_files,
                            fout=cst.neut_masks)

    # free up space from big arrays that are no longer needed
    del mcons, msubs, mcall

    # count substitutions at neutral sites to estimate neutral mutation rate
    subs = substitution_counts(nmask, 2e4, chrom)
    usegs, urates = prepare_umap(chrom, subs)

    # prepare bmaps
    bfs = cst.bmap_files.values()[0]
    seg, val = collect_bmaps(chrom, bfs)
    seg.append(usegs)
    val.append(urates)

    # combine maps along minimal segment intervals
    new_seg, new_val = combine_segments(chrom, seg, val)

    # get SNP table for neutral diversity calculations
    snp = np.column_stack(snpcount(cst.snp_files))

    # compress neutral data into segments
    comp = compress_data(nmask, snp, new_seg)

    # save segments as array
    np.savez_compressed(cst.seg_files, segs=new_seg)

    # save umap (last column of new_val)
    np.savez_compressed(cst.u_files, umap=new_val[:, -1])

    # save bmap array
    np.savez_compressed(cst.bs_files, bmap=new_val[:, :-1])

    # save compressed neutral data
    np.savez_compressed(cst.nt_files, neut=comp[:, :2])

    # save compressed divergence data
    np.savez_compressed(cst.ndiv_files, ndiv=comp[:, 2:])

    return None


def main_remote():
    if len(argv) != 3:
        print 'usage: pre_LH chrom init'
        exit(1)
    chrom = argv[1]
    init = argv[2]
    cst = ChromStruct(chrom=chrom, init=init)
    # cst = ChromStruct(chrom=chrom,
    #                   root=root_dir,
    #                   gmap='AA_Map_10kb',
    #                   refg='hg19',
    #                   mask='strictMask',
    #                   ncon='100vert',
    #                   cons='primate',
    #                   outg='rheMac3',
    #                   neut='YRI',
    #                   phase='phase3',
    #                   token='pr95.fixB',
    #                   focal='nonsyn',
    #                   rm_annos=default_filter_file,
    #                   bkgd_scale=100,
    #                   bs_annos=('primate_cons95_Segments',),
    #                   bmap_dir='bigB100',
    #                   cs_annos=(),
    #                   cmap_dir='cmaps',
    #                   bs_dfe=(default_fitness_effects,),
    #                   cs_dfe=(),
    #                   chroms=human_autosomes,
    #                   methods=('COBYLA', 'Powell'))

    # generate mask files
    mcons = conservation_mask(cst.wigz_ncons, chrom)
    msubs = substitution_mask(cst.axt_outgroup, chrom)
    mcall = sequence_mask(cst.call_mask, chrom)

    # create and save neutral mask
    nmask = neutrality_mask(mcons, mcall, msubs, filter_files=cst.filter_files,
                            fout=cst.neut_masks)

    # free up space from big arrays that are no longer needed
    del mcons, msubs, mcall

    # count substitutions at neutral sites to estimate neutral mutation rate
    subs = substitution_counts(nmask, 2e4, chrom)
    usegs, urates = prepare_umap(chrom, subs)

    # prepare bmaps
    bfs = cst.bmap_files.values()[0]
    seg, val = collect_bmaps(chrom, bfs)
    seg.append(usegs)
    val.append(urates)

    # combine maps along minimal segment intervals
    new_seg, new_val = combine_segments(chrom, seg, val)

    # get SNP table for neutral diversity calculations
    snp = np.column_stack(snpcount(cst.snp_files))

    # compress neutral data into segments
    comp = compress_data(nmask, snp, new_seg)

    # save segments as array
    np.savez_compressed(cst.seg_files, segs=new_seg)

    # save umap (last column of new_val)
    np.savez_compressed(cst.u_files, umap=new_val[:, -1])

    # save bmap array
    np.savez_compressed(cst.bs_files, bmap=new_val[:, :-1])

    # save compressed neutral data
    np.savez_compressed(cst.nt_files, neut=comp[:, :2])

    # save compressed divergence data
    np.savez_compressed(cst.ndiv_files, ndiv=comp[:, 2:])

    return None


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy'):
        main_local()
    else:
        main_remote()