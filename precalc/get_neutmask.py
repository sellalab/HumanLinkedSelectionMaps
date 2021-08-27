__author__ = 'davidmurphy'

from lh_inputs import ChromStruct, build_neutmask, argv


def main():
    if len(argv) != 3:
        print 'usage: get_neutmask <chrom> <init>'
        exit(1)
    cst = ChromStruct(chrom=argv[1], init=argv[2])
    build_neutmask(cst)


if __name__ == '__main__':
    main()
