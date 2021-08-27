import numpy as np
import sys

__author__ = 'davidmurphy'


def main():
    if len(sys.argv) != 2:
        print 'usage: txt2npy <txt_file>'
        exit(1)

    # remove suffix from txt_file and replace with .npy
    txt_file = sys.argv[1]
    npy_file = '.'.join(txt_file.split('.')[:-1]) + '.npy'

    # load data from txt_file and save a duplicate copy as in .npy format
    data = np.loadtxt(txt_file)
    np.save(npy_file, data)


if __name__ == '__main__':
    main()
