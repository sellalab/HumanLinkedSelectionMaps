import numpy as np
from itertools import izip

__author__ = 'davidmurphy'


class BkgdCalculator(object):
    """class that calculates b on sites x"""
    def __init__(self, coef, udel, anno, rcalc='haldane'):
        assert rcalc in ('linear', 'haldane')
        self.coef = coef  # selection coefficient
        self.udel = udel  # deleterious mutation rate
        self.anno = anno  # AnnoSegments instance
        self.rcalc = rcalc  # "haldane" or "linear" method to compute _rxy
        self.gcons = self.anno.rsites  # annotated sites in Morgans

    def __len__(self):
        """defined as the number of conserved sites effecting B"""
        return self.anno.num_bases

    def _rxy(self, gy, gx):
        """calculate linear distance or Haldane's r between gmap points"""
        assert self.rcalc in ('linear', 'haldane')
        m = abs(gy - gx)
        if self.rcalc == 'linear':
            r = m
        else:
            # use safe function expm1 to preserve precision (flip sign)
            r = -np.expm1(-2.0*m) / 2.0
        return r

    def bterm(self, r_xy):
        """return the little b term for distance r_xy"""
        return 1.0 / (self.coef * (1 + (self.rho * r_xy)) ** 2)

    def bkgd_sum(self, gy, gx, gpos=True, bsum=True):
        """exponential term for b at gy due to selection at points gx"""
        # convert to genetic map units if needed
        if not gpos:
            gy = self.gmap.interp_gpos(gy)
            gx = self.gmap.interp_gpos(gx)
        r_xy = self._rxy(gy, gx)
        if bsum:
            b = -self.udel * np.sum(self.bterm(r_xy))
        else:
            b = -self.udel * self.bterm(r_xy)
        return b

    def bkgd_sum_chrom(self, gy, gpos=True, bsum=True):
        """calculate B at y summing across conserved sites"""
        if not gpos:
            gy = self.gmap.interp_gpos(gy)
        gx = self.gcons  # use chromosome-wide annotation segments
        b = self.bkgd_sum(gy, gx, gpos=True, bsum=bsum)
        return b

    def dbterm(self, r_xy):
        """return the little b 1st derivative term for distance r_xy"""
        # return -2 * self.rho / (self.coef * (1 + (self.rho * r_xy)) ** 3)
        return 1.0 / (self.coef * (1 + (self.rho * r_xy)) ** 3)

    def bkgd_sum_der1(self, gy, gx, gpos=True, dbsum=True):
        """calculate the first derivative of b at gy"""
        # convert to genetic map units if needed
        if not gpos:
            gy, gx = self.gmap.interp_gpos(np.array([gy, gx]))
        r_xy = abs(gy - gx)  # distance between y and x in M
        # db = -self.udel * self.dbterm(r_xy)
        if dbsum:
            db = 2 * self.udel * self.rho * np.sum(self.dbterm(r_xy))
            # db = np.sum(db)
        else:
            db = 2 * self.udel * self.rho * self.dbterm(r_xy)
        return db

    def bkgd_sum_chrom_der1(self, gy, gpos=True, dbsum=True):
        """calculate B at y summing across conserved sites"""
        if not gpos:
            gy = self.gmap.interp_gpos(gy)
        gx = self.gcons  # use chromosome-wide annotation segments
        b = self.bkgd_sum_der1(gy, gx, gpos=True, dbsum=dbsum)
        return b

    def max_error_step(self, gy, gx, err):
        """maximum distance from gy where b changes by < err"""
        return err / self.bkgd_sum_der1(gy, gx)

    def bkgd_intgr(self, gy, gx1, gx2):
        """integral approximation to sum over points in conserved block"""
        dist = abs(gy-gx1), abs(gy-gx2)  # (safety measure)
        r1, r2 = min(dist), max(dist)
        trm1 = 1.0 / (self.rho*r1 + 1)
        trm2 = 1.0 / (self.rho*r2 + 1)
        b = -(trm1 - trm2) * self.udel / (1-self.coef)
        return b

    def bkgd_intgr_chrom(self, y, rates):
        """calculate B at y integrating across conserved segments"""
        # convert y to genetic map units
        gy = self.gmap.interp_gpos(y)
        rx = izip(self.anno.rcoords, rates)
        bterms = [self.bkgd_intgr(gy, ri, rj) * rt for ((ri, rj), rt) in rx]
        b = sum(bterms)
        # b = sum(self.bkgd_intgr(gy, ri, rj, rt) for ((ri, rj), rt) in rx)
        return b

    @property
    def gmap(self):
        return self.anno.gmap

    @property
    def chrom(self):
        return self.anno.chrom

    @property
    def rho(self):
        return (1 - self.coef) / self.coef
