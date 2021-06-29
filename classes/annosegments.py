import numpy as np
from geneticmap import GeneticMap

__author__ = 'davidmurphy'


class AnnoSegments(object):
    """a collection of segments for some annotation on a chromosome"""
    def __init__(self, gmap, afile=None, coords=None, one_based=False):
        assert isinstance(gmap, GeneticMap)
        self.gmap = gmap
        self.afile = afile
        self.coords = coords
        self._setcoords(one_based=one_based)

    def _setcoords(self, one_based):
        if self.afile:
            assert self.chrom in self.afile
            assert self.coords is None
            self.coords = np.loadtxt(self.afile, usecols=(1, 2))
        else:
            assert self.afile is None
        if one_based:
            self.coords[:, 0] -= 1  # convert from 1- to 0-based if flagged

    def nearest(self, ipos, start=True):
        imax = self.num_segs - 1
        if start:
            idx = np.minimum(imax, np.searchsorted(self.start, ipos))
        else:
            idx = np.minimum(imax, np.searchsorted(self.end, ipos))
        return idx

    def save_file(self, fout):
        slen = len(self.chrom)
        col1 = np.full(self.num_segs, self.chrom, dtype='S{}'.format(slen))
        arr = np.column_stack((col1, self.coords.astype(int).astype(str)))
        np.savetxt(fout, arr, fmt='%s\t%s\t%s')

    @property
    def chrom(self):
        return self.gmap.chrom

    @property
    def start(self):
        return self.coords[:, 0]

    @property
    def end(self):
        return self.coords[:, 1]

    @property
    def lengths(self):
        return self.end - self.start

    @property
    def sites(self):
        # site_array = np.zeros(shape=self.num_bases)
        # i = j = 0
        # for (start, end) in self.coords:
        #     j = i + (end - start)
        #     site_array[i:j] = np.arange(start, end)
        #     i = j
        return np.concatenate([np.arange(i, j) for (i, j) in self.coords])

    @property
    def rcoords(self):
        return self.gmap.interp_gpos(self.coords)

    @property
    def rstart(self):
        return self.gmap.interp_gpos(self.start)

    @property
    def rend(self):
        return self.gmap.interp_gpos(self.end)

    @property
    def rlengths(self):
        return self.rend - self.rstart

    @property
    def rsites(self):
        return self.gmap.interp_gpos(self.sites)

    @property
    def num_segs(self):
        return len(self.coords)

    @property
    def num_bases(self):
        return self.lengths.sum()
