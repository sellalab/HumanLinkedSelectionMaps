__author__ = 'davidmurphy'


from sys import argv
from classes.runstruct import chromosome_length, np


def main():
    if len(argv) != 2:
        print 'usage: fix_bmaps <f_bmap>'
        exit(1)

    # get the bmap file from command line
    f_bmap = argv[1]

    # get chrom from bmap file name
    chrom = f_bmap.split('/')[-1].split('.')[1]

    # load file with header and data in separate lists
    head = []
    data = []
    with open(f_bmap, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # save header strings
                head.append(line)
            else:
                # convert row data into integers
                row_data = map(int, line.split())
                data.append(row_data)

    # convert data into an array
    data = np.array(data, dtype=int)

    # calculate length and compare to real length
    if data[:,1].sum() != chromosome_length(chrom):
        # add additional site at the end if length mismatch
        data[-1,1] += 1
    # single base pair should fix the gaps
    assert data[:,1].sum() == chromosome_length(chrom)

    # re-write the file
    with open(f_bmap, 'w') as f:
        for line in head:
            f.write(line)
        for (b, l) in data:
            row = '{} {}\n'.format(b, l)
            f.write(row)


if __name__ == '__main__':
    main()
