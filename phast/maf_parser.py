#!/usr/bin/python

from classes.runstruct import ChromStruct
from maf_tools import Parser, Splitter, PhastCaller
from sys import argv, stdin
from gzip import open as zopen
import os
import re

__author__ = 'davidmurphy'

# path for downloading 100 vertebrate MAF files from UCSC:
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/maf/chr2.maf.gz


def main_local():
    cst = ChromStruct(chrom='chr21')
    odir = cst.data + '/phast'
    sp_str = 'euarchontoglires:hg19,panTro4,gorGor3,ponAbe2,nomLeu3,rheMac3,macFas5,papHam1,chlSab1,calJac3,saiBol1,otoGar3,tupChi1,speTri2,jacJac1,micOch1,criGri1,mesAur1,mm10,rn5,hetGla2,cavPor3,chiLan1,octDeg1,oryCun2,ochPri3'
    tax, spc = sp_str.split(':')
    splitter = Splitter('chr21', tax, spc, odir)
    # split until the entire MAF file is complete
    while True:
        try:
            splitter.split_file().next()
        except StopIteration:
            break

    # # chromstruct
    # cst = ChromStruct(chrom='chr1')
    #
    # # taxa re pattern
    # txa = re.compile('[a-zA-Z]+\d+')
    #
    # # open phylo tree file and extra taxa list to make species string
    # sfile = '{}/data/phast/46way.nh'.format(cst.root)
    # with open(sfile, 'r') as f:
    #     txt = f.read()
    #     spc = txa.findall(txt)
    #     # trim species to placental mammals only
    #     spc = spc[:35]
    # species = ','.join(spc)
    #
    # # create a Parser to process lines of maf file into blocks
    # parser = Parser(species_string=species)
    #
    # # # test maf file
    # # maf_file = '{}/data/phast/chr22.maf.gz'.format(cst.root)
    # #
    # # # counter for syntenic bases
    # # n_syntenic = 0
    #
    # # open a file to write blocks to
    # with open('46way.syntenic.blocks.txt', 'w') as f:
    #     while True:
    #         try:
    #             # create blocks from stdin lines of MAF data (from pipe)
    #             blk = parser.get_block(stdin, as_string=False).next()
    #             # if block has data from every species, take hg19 coords
    #             if len(blk) > len(spc):
    #                 seq = blk[1]
    #                 assert seq.species == 'hg19'
    #                 f.write('{} {} {}\n'.format(seq.chrom, seq.start, seq.size))
    #         except StopIteration:
    #             break

    # model = out_dir + '/primate.noncons.mod'
    # caller = PhastCaller(chrom, model, 1e5, init_string, maf_file)
    # caller.send_jobs(call_limit=1)

    # chrom = 'chr1'
    # taxon = 'primate'
    # species = 'hg19,panTro4,gorGor3,ponAbe2,nomLeu3,rheMac3,macFas5,papHam1,chlSab1'
    # splitter = Splitter(chrom, taxon, species, 1e5, out_dir, maf_file)
    # n = 0
    # blocks = []
    # while True:
    #     try:
    #         print splitter.split_file().next()
    #     except StopIteration:
    #         break


def main_remote():
    if len(argv) != 3:
        print 'usage: maf_parser <chrom> <init>'
        exit(1)

    # set output file path
    odir = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/mafs'
    # get chrom, init file from command line
    ch, init = argv[1:3]
    with open(init, 'r') as f:
        tax, spc = f.read().split(':')

    # initialize the splitter. split until the entire MAF file is complete
    splitter = Splitter(ch, tax, spc, odir)
    while True:
        try:
            splitter.split_file().next()
        except StopIteration:
            break


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
