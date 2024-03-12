from classes.annosegments import AnnoSegments
from classes.runstruct import ChromStruct
from classes.geneticmap import GeneticMap
from precalc.bkgd_grid import bkgd_range
from data_tools import calculate_error
from itertools import izip
import numpy as np
import os

__author__ = 'davidmurphy'


def main():
    cst = ChromStruct('chr1')
    an = cst.bs_annos[0]
    num_divisions = 4
    # for ch in cst.chroms:
    for ch in ['chr22']:
        # initalize new map and annotation classes per chrom
        cst.chrom = ch
        gmp = GeneticMap(ch, cst.gmap_files, 0.01)
        ano = AnnoSegments(ch, gmp, cst.bs_targets[an])

        # find M/bp for each segment
        dp_dm = ano.rlengths / ano.lengths
        med_rate = np.median(dp_dm)

        # reorganize coords by splitting those that have low M/bp rates
        new_coords = []
        for (coord, rate) in izip(ano.coords.astype(int), dp_dm):
            start, end = coord
            num_bases = end - start
            if rate < med_rate:
                # limit new segment lengths at 1bp
                new_len = max(1, int(num_bases / num_divisions))
                # determine how many segments are needed for subdivision
                num_segs = int(np.ceil(float(num_bases) / float(new_len)))
                # add each new segment sequentially
                for _ in range(num_segs-1):
                    s_i = start
                    s_j = start + new_len
                    new_coords.append((s_i, s_j))
                    start = s_j

                # add remaining segment if original was irregular length
                if start != end:
                    new_coords.append((start, end))

            else:
                # otherwise leave the segment as is
                new_coords.append((start, end))

        # check that the total number of bases is the same
        new_coords = np.array(new_coords, dtype=int)
        assert ano.num_bases == np.sum(new_coords[:, 1]-new_coords[:, 0])

        # save the new split annotations in a modified directory
        new_afile = ano.afile.replace('.bed', '_split4.bed')
        slen = 'S{}'.format(len(ano.chrom))
        col_1 = np.full(len(new_coords), ano.chrom, dtype=slen)
        aout = np.column_stack((col_1, new_coords.astype(str)))
        np.savetxt(new_afile, aout, fmt='%s\t%s\t%s')


def profile_map_errors(bdir, eps):
    """get information from regions with large "exact" McVicker B errors"""
    anno = 'primate_cons95_Segments'
    token = 'pr95.cleanrun'
    tval = 10**-4.5
    suff = 'edge'

    cst = ChromStruct('chr1', bdir=bdir, tkn=token)
    for ch in cst.chroms:
        # update chromstruct to update file paths
        cst.chrom = ch
        bfile = cst.bxct_file(anno, tval, suff)
        bx, b_map, b_xct = np.load(bfile).T

        # select points where abs(rel_err) exceeds tolerance set by epsilon
        rel_err = calculate_error(b_xct, b_map, rel=True, fabs=True)
        r_pts = bx[(rel_err > 2*eps)]

        # analyze data within the effective BS range surrounding error points
        gmp = GeneticMap(ch, cst.gmap_files)
        ano = AnnoSegments(gmp, afile=cst.bs_target(anno))
        idx = np.zeros(shape=(len(r_pts), 2))
        for (i, r) in enumerate(r_pts):
            si, sj = bkgd_range(r, cst, ano, tval, eps)
            pi, pj = ano.sites[si], ano.sites[sj]
            idx[i, :] = pi, pj

        return idx


def remove_smallest(pct):
    cst = ChromStruct('chr1')
    an = 'primate_cons95_Segments'
    new_an = 'primate_cons95_Segments_trim{}'.format(pct)
    new_dir = '{}/data/bsanno/{}'.format(cst.root, new_an)
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)

    for ch in cst.chroms:
        cst.chrom = ch
        gmp = GeneticMap(ch, cst.gmap_files)
        ano = AnnoSegments(gmp, cst.bs_target(an))
        cutoff = np.percentile(ano.rlengths, pct)
        msk = (ano.rlengths > cutoff)
        print float(msk.sum()) / float(ano.num_segs)
        ano.coords = ano.coords[msk]
        fout = ano.afile.replace(an, new_an)
        ano.save_file(fout)


if __name__ == '__main__':
    remove_smallest(10)
