__author__ = 'davidmurphy'


from sys import argv
from classes.runstruct import chromosome_length


def fill_file_gap(f_name):
    # AA_Map_10kb.chr10.h1hesc_enh.t0.00003162.bkgd
    # get chrom from filename
    ch = f_name.split('/')[-1].split('.')[1]

    # read file contents into 3 lists
    head = []
    bval = []
    lens = []
    with open(f_name, 'r') as f:
        for line in f:
            # store header separately
            if line.startswith('#'):
                head.append(line)
            else:
                b, l = line.split()
                bval.append(b)
                lens.append(int(l))

    # calculate sum of segments
    seglen = sum(lens)

    # check for a gap
    gap = chromosome_length(ch) - seglen

    # if there is a gap, fill it in the final segment and rewrite the file
    if gap > 0:
        print '{} gap size: {}'.format(f_name, gap)
        lens[-1] += gap
        f_save = f_name  # + '.fixed'
        with open(f_save, 'w') as f:
            for line in head:
                f.write(line)
            for (b, l) in zip(bval, lens):
                line = '{} {}\n'.format(b, l)
                f.write(line)

    return None


def main():
    if len(argv) != 2:
        print 'usage: fill_tail_gaps <filename>'
        exit(1)
    fill_file_gap(argv[1])


if __name__ == '__main__':
    main()
