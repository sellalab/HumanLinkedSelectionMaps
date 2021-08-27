import os
import sys
import gzip


__author__ = "davidmurphy"


def main():
    """
    collect each small wig file for a given chrom-species pair
    into one zipped file and then erase the small files
    """
    
    if len(sys.argv) != 3:
        print "usage: <chrom> <tree>"
        exit(1)

    # load specific chromosome and species tree
    chrom = sys.argv[1]
    tree  = sys.argv[2]
    
    # master file to write to compressed
    out_file = "zips/{tr}.{ch}.wig.gz".format(tr=tree, ch=chrom)

    # clear pre-existing file if it exists so that we don't double append anything
    if os.path.isfile(out_file):
        os.remove(out_file)
    writable = gzip.open(out_file, "a")

    # loop through all small "n_<X>" prefixed wig files, write them to a single file, then erase them
    n = 1
    cur_file = "wigs/n_{num}.{tr}.{ch}.wig".format(num=n, tr=tree, ch=chrom)
    while os.path.isfile(cur_file):
        # print cur_file
        with open(cur_file, "r") as f:
            data = f.read()
        writable.write(data)
        os.remove(cur_file)
        n += 1
        cur_file = "wigs/n_{num}.{tr}.{ch}.wig".format(num=n, tr=tree, ch=chrom)

    writable.close()

    return None


if __name__ == "__main__":
    main()
