#!/ifs/data/c2b2/gs_lab/shared/software_hpc/anaconda/bin/python


__author__ = 'davidmurphy'


import os
import subprocess
from sys import argv
from classes.runstruct import ChromStruct, np


def chlen_dict():
    """get a chrom-length dictionary"""
    # path to hg19 chrom length file
    if os.getcwd().startswith('/Users/davidmurphy'):
        f_chlen = '/Users/davidmurphy/GoogleDrive/linked_selection/data' \
                  '/ch_features/hg19_chroms.txt'
    else:
        f_chlen = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data' \
                  '/coords/chr_len_all.txt'
    # open file and convert chrom, length pairs into a dictionary
    ch_dict = {}
    with open(f_chlen, 'r') as f:
        f.readline()  # skip comment line
        for line in f:
            chrom, nbase = line.split()
            ch_dict[chrom] = int(nbase)

    return ch_dict


def main_0():
    """call ffit_windows over length of given chromosome"""
    if len(argv) != 4:
        print 'usage: call_ffit <ch> <size> <win>'
        exit(1)

    # get chrom, total segment size and window size from command line
    ch = argv[1]
    size, win = map(int, map(float, argv[2:]))

    # get chrom length
    clen = chlen_dict()[ch]

    # use primate as hardcoded default multiple alignment
    ma = 'primate'

    # format the shell command
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    prog = '{}/ffit_windows.sh'.format(phast_dir)
    cfmt = '{pr} {ch} {ma} {st} {en} {wn}'

    # iterate over chunks of the chromosome and subdivide them into 20kb
    for start in range(0, clen, int(size)):
        # format and execute the shell command on each chunk
        end = start+int(size)
        cmd = cfmt.format(pr=prog, ch=ch, ma=ma, st=start, en=end, wn=win)
        subprocess.call(cmd, shell=True)


def main_1():
    """call ffit_windows over length of given chromosome"""
    if len(argv) != 4:
        print 'usage: call_ffit <ch> <size> <win>'
        exit(1)

    # fix the species number in the alignment to be used
    n_spc = 8

    # get chrom, total segment size and window size from command line
    ch = argv[1]
    # nspec, size, win = map(int, map(float, argv[2:]))
    size, win = map(int, map(float, argv[2:]))

    # get compressed neutral aligned FASTA length from neutmask
    cst = ChromStruct(chrom=ch)
    new_tkn = '.aligned.{}.nmsk'.format(n_spc)
    f_mask = cst.neut_masks.replace('.nmsk', new_tkn)
    nmsk = np.load(f_mask)['neutmask']
    clen = nmsk.sum()

    # format the shell command
    phast_dir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast'
    # prog = '{}/ffit_windows.sh'.format(phast_dir)
    # cfmt = '{pr} {ch} {ns} {st} {en} {wn}'
    prog = '{}/sh/cphyloseg.sh'.format(phast_dir)
    cfmt = '{pr} {ch} {st} {en} {wn}'

    # iterate over chunks of the chromosome and subdivide them into 20kb
    for start in range(0, clen, int(size)):
        # format and execute the shell command on each chunk
        end = min(start+int(size), clen)
        cmd = cfmt.format(pr=prog, ch=ch, st=start, en=end, wn=win)
        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main_1()
