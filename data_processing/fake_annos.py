from data_tools import randint_unique
from classes.annosegments import AnnoSegments
from classes.runstruct import ChromStruct
from classes.geneticmap import GeneticMap
from itertools import izip
import numpy as np
import os

__author__ = 'davidmurphy'


def generate_segments(rmin, rmax, gmap, n_bp):
    """create random segments in regions with M/bp between rmin and rmax"""
    blocks = np.where((gmap.rate >= rmin) & (gmap.rate <= rmax))[0]
    blocks = blocks[blocks < (len(gmap.pos) - 2)]
    ranges = [(gmap.pos[i], gmap.pos[i+1]) for i in blocks]
    pct = 1.0 * n_bp / sum(b-a for (a, b) in ranges)

    fake_segments = []
    for r in ranges:
        total_segs = int((r[1]-r[0]-10) / 10)
        size = int(pct * total_segs)
        idx = randint_unique(size, total_segs)
        fake = [(r[0] + i*10, r[0] + i*10 + 10) for i in idx]
        if any(fake):
            fake_segments.append(fake)

    segments = np.concatenate(fake_segments)
    return segments


def main():
    an = 'primate_cons95_Segments'
    # new_an = 'low_rec_rate'
    # new_an = 'high_rec_rate'
    new_an = 'rm_norec'
    cst = ChromStruct('chr1')

    # for ch in cst.chroms:
    for ch in ['chr22']:
        cst.chrom = ch
        gmp = GeneticMap(ch, cst.gmap_files, 0.01)
        ano = AnnoSegments(gmp, afile=cst.bs_target(an))
        r = np.searchsorted(gmp.pos, ano.start) == \
            np.searchsorted(gmp.pos, ano.end)

        ano.coords = ano.coords[ano.rstart != ano.rend]

        # find median rate in chrom map and take scores above or below
        # med = np.median(gmp.rate)
        # low
        # coords = generate_segments(0.0, med, gmp, ano.num_bases)
        # high
        # coords = generate_segments(med, gmp.rate.max(), gmp, ano.num_bases)

        new_dir = '{}/data/bsanno/{}'.format(cst.root, new_an)
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)
        fout = ano.afile.replace(an, new_an)
        ano.save_file(fout)


if __name__ == '__main__':
    main()