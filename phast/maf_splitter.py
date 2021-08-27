from maf_tools import Parser, Splitter
from sys import argv
import os

__author__ = "davidmurphy"


def main_remote():

    if len(argv) != 5:
        print "usage: maf_splitter <maf_file> <config_file> <model> <out_pref>"
        exit(1)

    # use a single config file to get a list of species subset files and the chromosome name for the new MAFs
    maf_file    = argv[1]
    config_file = argv[2]
    model_file  = argv[3]
    out_pref    = argv[4]

    # make sure the chromosome matches the filename
    assert out_pref in maf_file

    # maf_file = "/Users/davidmurphy/GoogleDrive/linked_selection/data/conservation/sample_maf.maf"
    # config_file = "/Users/davidmurphy/GoogleDrive/linked_selection/data/conservation/species_cfg.txt"
    # model_file = "/Users/davidmurphy/GoogleDrive/linked_selection/data/conservation/hg19.100way.phastCons.mod"
    # out_pref = "chr1"


def main_local():
    pass


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
