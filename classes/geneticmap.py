import os
import numpy as np
from classes.runstruct import chromosome_length, root_dir


__author__ = 'davidmurphy'


class GeneticMap(object):
    """a class to represent genetic maps and perform some basic functions"""
    def __init__(self, chrom, gfile, gscale=1.0):
        assert chrom in gfile
        self.chrom = chrom
        self.gfile = gfile
        self.gscale = gscale
        self._gmap = self._getmap()

    def __len__(self):
        return len(self.pos)

    def _getmap(self):
        """load genetic map, format=[position, cM/Mb, cM] using gscale"""
        assert isinstance(self.gscale, float)
        gmap = np.loadtxt(self.gfile)
        if self.gscale != 1.0:
            gmap[:, 1:] *= self.gscale
        return gmap

    def reload_map(self):
        """reset map values from map file"""
        self._gmap = self._getmap()

    def gmap_mask(self, ipos):
        """remove ipos outside genetic map pos limits"""
        gmask = (ipos > self.pos[0]) & (ipos < self.pos[-1])
        return ipos[gmask]

    def interp_gpos(self, ipos):
        """convert physical map position to genetic map position"""
        return np.interp(ipos, self.pos, self.gpos)

    def interp_rate(self, ipos):
        """return the recombination rate at gmap bin corresponding to pos"""
        idx = np.minimum(len(self)-1, np.searchsorted(self.pos, ipos)) - 1
        return self.rate[idx]

    def interp_pos(self, igpos):
        """approximate interpolation from genetic position to physical"""
        return np.interp(igpos, self.gpos, self.pos).astype(int)

    def close_gaps(self):
        """stretch gmap to tile the entire reference chromosome"""
        pos, rate, gpos = self.gmap.T
        # start map from position 1
        if pos[0] > 1:
            # stretch map down to 1
            pos = np.concatenate(([1], pos))
            # assign first genetic position to be 0.0
            gpos = np.concatenate(([0.0], gpos))
            # calculate the recombination rate from the new first position
            d_r = (1e6 * gpos[1]) / pos[1]
            rate = np.concatenate(([d_r], rate))
        # end map at position chromosome_length
        if pos[-1] < self.chlen:
            # stretch map up to chromosome length
            pos = np.concatenate((pos, [self.chlen]))
            # assign the same genetic position as the last gpos
            gpos = np.concatenate((gpos, [gpos[-1]]))
            # no recombination rate since gpos did not chance
            rate = np.concatenate((rate, [0]))
        assert np.sum(pos[1:] - pos[:-1]) + 1 == self.chlen
        self._gmap = np.column_stack((pos, rate, gpos))

    def save_map(self, new_fname=None):
        """save the current working copy of the map to map file"""
        # if no new map name is provided, add "COPY" to original
        if new_fname is None:
            fout = self.gfile.replace('.txt', '.COPY.txt')
        else:
            fout = new_fname
        # format output and map header
        fmt = '%d\t%.14g\t%.14g'
        header = 'position\trec_rate(cM/Mb)\tgenetic_map(cM)'
        np.savetxt(fout, self.gmap, fmt=fmt, header=header)

    @property
    def gmap(self):
        return self._gmap

    @property
    def pos(self):
        return self.gmap[:, 0]

    @property
    def rate(self):
        return self.gmap[:, 1]

    @property
    def gpos(self):
        return self.gmap[:, 2]

    @property
    def gmap_limits(self):
        return self.gpos.min(), self.gpos.max()

    @property
    def chlen(self):
        return chromosome_length(self.chrom)

    @property
    def viable_5prime(self):
        """return the 5' position where the map is putatively viable"""
        imax = np.where(self.gpos == self.gpos[self.gpos > 0].min())[0].max()
        return self.pos[imax].astype(int)

    @property
    def viable_3prime(self):
        """return the 3' position where the map is putatively viable"""
        jmin = np.where(self.gpos == self.gmap_limits[1])[0].min()
        return self.pos[jmin].astype(int)


def rescale_map(gmap, scale):
    """rescale genetic distances"""
    assert isinstance(scale, float)
    if scale != 1.0:
        gmap[:, 1:] *= scale
    return gmap


def anjali2gmap(chrom, anjali_file, map_id, res):
    """
    Create a new
    :param chrom: map chromosome
    :param anjali_file: anjali map file to open
    :param map_id: name of the map column to take (must be in anjali file)
    :param res: resolution in bp for sampling recombination rates
    """
    # create path to new output file, create map directory if needed
    fdir = '{}/data/maps/{}'.format(root_dir, map_id)
    fout = '{}/{}_{}.txt'.format(fdir, chrom, map_id)
    if not os.path.isdir(fdir):
        os.mkdir(fdir)

    # write the header to the outfile
    fh_out = open(fout, 'w')
    fh_out.write('#{} taken from Anjali hg19 liftover\n'.format(map_id))
    fh_out.write('# position cM/Mb cM\n')

    # set initial position to 1, genetic map position to 0.0
    cur_pos = new_pos = 1
    cur_gpos = new_gpos = 0.0
    cm_mb = 0

    # # record initial pos, cmmb and gpos
    # rec = '{} {:.14g} {:.14g}\n'.format(cur_pos, 0.0, cur_gpos)
    # fh_out.write(rec)

    # open anjali map file
    with open(anjali_file, 'r') as f:
        # get map id index from the header before reading the map
        header = f.readline().replace('"', '').split()
        idx = header.index(map_id)

        # iterate over the map and take genetic map positions from idx column
        for line in f:
            line = line.split()
            new_pos = int(line[0])
            new_gpos = float(line[idx])
            d_pos = new_pos - cur_pos  # change since last position

            # create new record when delta_pos > res
            if d_pos >= res:
                d_mb = 1e-6 * d_pos  # delta_pos in Mb
                d_cm = new_gpos - cur_gpos
                cm_mb = d_cm / d_mb  # local recombination rate in cM/Mb
                rec = '{} {:.14g} {:.14g}\n'.format(cur_pos, cm_mb, cur_gpos)
                fh_out.write(rec)

                # reset current position and genetic position
                cur_pos = new_pos
                cur_gpos = new_gpos

    # create a final record to complete the chromosome
    if new_pos != cur_pos:
        d_pos = new_pos - cur_pos  # change since last position
        d_mb = 1e-6 * d_pos  # delta_pos in Mb
        d_cm = new_gpos - cur_gpos
        cm_mb = d_cm / d_mb  # local recombination rate in cM/Mb
    last_pos = chromosome_length(chrom)  # set new rate up to the final pos
    rec = '{} {:.14g} {:.14g}\n'.format(last_pos, cm_mb, new_gpos)
    fh_out.write(rec)
    fh_out.close()

    return None
