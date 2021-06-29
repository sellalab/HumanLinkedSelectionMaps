import numpy as np
from geneticmap import GeneticMap

__author__ = 'davidmurphy'


class AnnoPoints(object):
    """a collection of chromosomal positions for an annotation"""
    def __init__(self, gmap, afile):
        assert isinstance(gmap, GeneticMap)
        self.gmap = gmap
        self.afile = afile
        if afile.endswith('.npz'):
            self.points = np.load(afile)['pos']
        else:
            self.points = np.load(afile)

    def __len__(self):
        return len(self.points)

    @property
    def chrom(self):
        return self.gmap.chrom

    @property
    def rpoints(self):
        """return points in genetic map units"""
        return self.gmap.interp_gpos(self.points)
