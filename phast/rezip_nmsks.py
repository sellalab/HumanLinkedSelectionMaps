__author__ = 'davidmurphy'


import numpy as np
from sys import argv


def rezip(f):
    a = np.load(f)['neutmask']
    np.savez_compressed(f, neutmask=a)


def main():
    if len(argv) != 2:
        print 'usage: rezip_nmsks <file>'
        exit(1)
    rezip(argv[1])


if __name__ == '__main__':
    main()
