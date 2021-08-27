__author__ = 'davidmurphy'

import numpy as np
from sys import argv
from gzip import open as zopen
from classes.runstruct import ChromStruct, root_dir
from data_processing.data_tools import randint_unique


def build_snpcounts(cst, sim_file):
    """
    Using simulated polymorphism data with sample size = 2, which has been
    summaraized into genomic segments (counts of hom, het pairs), create
    a frq.count format file type, wich bases randomly assigned (they will not
    be used for the analyses).
    :type cst: ChromStruct
    :param sim_file: path to simulation file used to build SNP data
    """
    # check that cst and sim_file match
    assert cst.chrom in sim_file

    # load neutral mask file
    nm = np.load(cst.neut_masks)['neutmask']
    # load neutral polymorphism from the simulation file
    nt = np.load(sim_file)['nt'].astype(int)
    # check that neutral sites sum to total hom+het pairs
    print 'nt={} nm>0={}'.format(np.sum(nt), np.sum(nm>0))
    # assert np.sum(nt) == np.sum(nm > 0)

    # get the positions of neutral sites
    neut_pos = np.where(nm > 0)[0] + 1
    # load neutral site count per segment of chrom
    ns = np.load(cst.dv_files)['dv'][:, 0].astype(int)
    # check the sites marked neutral in mask sum to neutral sites/segment
    print 'ns={} neut_pos={}'.format(np.sum(ns), len(neut_pos))
    assert np.sum(ns) == len(neut_pos)

    # rename snp file
    snp_file = cst.snp_files.replace('.phase3', '.sim.rand00.n2.ds00')

    # create template for SNP data lines
    # 1       10642   2       216     G:214   A:2
    snpstr = '{}\t{}\t2\t2\tA:{}\tG:{}\n'

    # write new data to file
    with zopen(snp_file, 'w') as f:
        # use snp site count per segment for POSITION index
        i = 0
        # use enumerate to get SEGMENT index
        for (idx, n) in enumerate(ns):
            # skip segments with no neutral sites
            if n == 0:
                continue

            # adjust the index
            j = i + n

            # get HET pairs for current segment
            nhet = nt[idx, 1]
            if nhet > 0:
                # check that nhet <= n
                assert nhet <= n
                # create an empty array for simulated het sites
                data = np.zeros(shape=(nhet, 3), dtype=int)
                # select random sites indices for the current segment
                ridx = randint_unique(nhet, n)
                # get the SNP positions at selected sites
                data[:, 0] = neut_pos[i:j][ridx]
                # set ref/alt to 1/1 at selected sites
                data[:, 1:] = [1, 1]

                # write simulated data lines to new file
                for (pos, ref, alt) in data:
                    ch = cst.chrom[3:]
                    line = snpstr.format(ch, pos, ref, alt)
                    f.write(line)

            # update i value
            i = j


def main():
    if len(argv) != 4:
        print 'usage: simulate_snpcounts <chrom> <init> <sim_file>'
        exit(1)
    else:
        chrom, init, sim_file = argv[1:]
        cst = ChromStruct(chrom=chrom, init=init)
        build_snpcounts(cst=cst, sim_file=sim_file)


if __name__ == '__main__':
    main()
