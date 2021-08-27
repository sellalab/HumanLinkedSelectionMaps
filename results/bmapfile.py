#!/usr/bin/python
from classes.runstruct import RunStruct, np, os

__author__ = 'davidmurphy'


def write_bmapfile(rst, chrom='', path=''):
    """
    Write a bmap file using McVicker formatting for chrom if specified otherwise all 22 autosomes. Store the files in a
    set path if specified, otherwise generate a default path.
    :param rst: RunStruct from a final inference run
    :param chrom: chromosome like 'chr1'
    :param path: optional directory to write the bmap files to
    :type rst: RunStruct
    :type chrom: str
    :type path: str
    """
    # check the write-path for the bmaps
    if not path:
        # default path name generated
        path = '{}/result/bmaps/{}_{}'.format(rst.root, rst.token, rst.vars.timestamp)
    if not os.path.isdir(path):
        # path written if it doesn't yet exist
        os.mkdir(path)

    # write a detailed header for the map
    header = ''

    # create a bmap array and write the file for "chrom" or all chromosomes
    if chrom:
        chroms = [chrom]
    else:
        chroms = rst.chroms
    t = []
    for ch in chroms:
        filename = '{}/{}.bkgd'.format(path, ch)
        arr = rst.get_bmap(chrom=ch)
        # arr[:, 0] *= 10  # convert to 0-1000 scale following McVicker formatting
        arr = arr.astype(int)  # convert to int
        # collapse segments where the values are unchanged in the map
        # collapsed = [list(arr[0])]  # start with the initial b value and segment length
        # for (b, seg) in arr[1:]:
        #     if b == collapsed[-1][0]:
        #         # the b value has not changed, just increment the segment
        #         collapsed[-1][1] += seg
        #     else:
        #         # new b value encountered, start a new segment
        #         collapsed.append([b, seg])
        # # save the collapsed bmap to a file
        # with open(filename, 'w') as f:
        #     f.write(header)
        #     f.write('\n'.join('{:d} {:d}'.format(b, seg) for (b, seg) in collapsed))
        t.append(arr)
    t = np.concatenate(t)
    import matplotlib.pyplot as plt
    plt.figure(figsize=(16, 9))
    plt.hist(x=t[:, 0], bins=t[:,0].max() - t[:,0].min(), weights=t[:,1])
    plt.show()
        # plt.plot(np.cumsum(arr[:,1]), arr[:,0])
        # x, y = np.array(collapsed, dtype=int).T
        # plt.plot(np.cumsum(y), x)
        # plt.show()

    return None


def main():
    f = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/yri-div/' \
    '08-prim-95-cons-div_BS1_CS0_161206111835_final.txt'
    r = RunStruct(init=f)
    write_bmapfile(rst=r)

    # # parse args from command line
    # parser = ArgumentParser(description='Call the full inference pipeline or individual functions')
    # parser.add_argument('init_file', type=str, help='path to a RunStruct init file')
    # parser.add_argument('-chrom', type=str, default='', help='optional chromosome name i.e."chr1"')
    # parser.add_argument('-path', type=str, default='', help='optional path to save bmaps to')
    # args = parser.parse_args()
    # # init the rst
    # rst = RunStruct(init=args.init_file)
    # write_bmapfile(rst=rst, chrom=args.chrom, path=args.path)
    return None

if __name__ == '__main__':
    main()
