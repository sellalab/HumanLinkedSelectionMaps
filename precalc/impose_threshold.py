__author__ = 'davidmurphy'


from sys import argv
from precalc.fill_tail_gaps import fill_file_gap


def impose_bmap_threshold(f_name, threshold=500):
    """join all segments above a given threshold to and give threshold value"""

    # read file contents into 3 lists
    head = []
    bval = []
    segs = []
    with open(f_name, 'r') as f:
        for line in f:
            # store header separately
            if line.startswith('#'):
                head.append(line)
            else:
                b, s = map(int, line.split())
                # apply threshold to b
                b = min(threshold, b)
                # if previous b value is the same after thresholding, join segs
                if len(bval) and bval[-1] == b:
                    segs[-1] += s
                else:
                    bval.append(int(b))
                    segs.append(int(s))
    # rewrite file
    f_save = f_name #+ '.thresh'
    with open(f_save, 'w') as f:
        for line in head:
            f.write(line)
        for (b, s) in zip(bval, segs):
            line = '{} {}\n'.format(b, s)
            f.write(line)

    return None





# f = '/Users/davidmurphy/GoogleDrive/linked_selection/precalc/chromHMM/AA_Map_10kb.chr1.h1hesc_txn.t0.00003162.bkgd'
# impose_bmap_threshold(f)


def main():
    if len(argv) != 2:
        print 'usage: impose_threshold <filename>'
        exit(1)
    impose_bmap_threshold(argv[1])
    # fill_file_gap(argv[1])


if __name__ == '__main__':
    main()
