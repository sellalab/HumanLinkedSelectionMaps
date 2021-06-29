#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 04:34:29 2018

@author: davidmurphy
"""

from collections import OrderedDict
from itertools import izip
import numpy as np
import re


class BkgdMap(object):
    """a class to represent bmaps and perform basic functions"""

    def __init__(self, bkgd, segs):
        """NOTE: BkgdMap always takes little bkgd values"""
        self._bkgd = bkgd
        self._segs = segs
        
    def interp_idx(self, pos):
        """get the map index for a given position"""
        imax = self.segs.size - 1
        idx = np.minimum(np.searchsorted(self.pnts, pos), imax)
        return idx

    def interp_bkgd(self, pos):
        """get b at given position on the map"""
        idx = self.interp_idx(pos)
        return self.bkgd[idx]

    def write_bmap(self, f_out, bscale, info=None, expb=False):
        """save bmap in a standardized text file format"""
        # discretize b values into positive integer bins
        b = self.ebkgd if expb else -self.bkgd
        b_int = np.floor(b * bscale + 0.5)
        # if info dict is provided, set BKGD_SCALE to bscale
        if info is not None:
            info['BKGD_SCALE'] = '{:.6f}'.format(bscale)
        else:
            info = dict(BKGD_SCALE='{:.6f}'.format(bscale))
        # write a file header followed by the formatted bmap data
        with open(f_out, 'w') as f:
            for (key, value) in info.items():
                f.write('#{}={}\n'.format(key, value))
            for (b, s) in izip(b_int, self.segs):
                f.write('{:.0f} {:.0f}\n'.format(b, s))

    @property
    def bkgd(self):
        return self._bkgd
    
    @property
    def ebkgd(self):
        return np.exp(self.bkgd)
   
    @property
    def segs(self):
        return self._segs

    @property
    def pnts(self):
        return np.cumsum(self.segs)
 
    
class BkgdMapReader(object):
    """a class to represent bmaps and perform basic functions"""
    # TODO: implement flexible map reading for big B scaled maps
    def __init__(self, bfile=None, safemode=True):
        self._info = OrderedDict()
        self._data = []
        if bfile:
            self.read_bfile(bfile, safemode)
        
    def read_bfile(self, bfile, safemode=True):
        """process bmap from file. call iteratively to merge partitions"""
        # last partition index (starts at -1 before there is any info)
        last_pidx = self.info['CHR_PARTITION_INDEX'] if self.info else -1
        with open(bfile, 'r') as f:
            # tally the total segment length for each read file
            s_sum = 0
            for line in f:
                # process header
                if line.startswith('#'):
                    self._hdrinfo(line)
                # process data: merge adjacent segs with matching b values
                else:
                    l = line.split()
                    if len(l) != 2:
                        # msg = 'pidx={}: incomplete line encountered'
                        msg = 'file={} incomplete line encountered'
                        curr_pidx = self.info['CHR_PARTITION_INDEX']
                        raise ValueError(msg.format(bfile))
                    try:
                        b, s = map(float, l)
                    except ValueError:
                        msg = 'file={} corrupted data encountered'
                        raise IOError(msg.format(bfile))
                    if (self.data) and (self.data[-1][0] == b):
                        self._data[-1][1] += s
                    else:
                        self._data.append([b, s])
                    s_sum += s

        # PATCH: THIS LINE FIXES SMALL GAPS AT THE END OF BMAPS (up to 2bp)
        try:
            p_len = self.info['CHR_PARTITION_LENGTH']
        except KeyError:
            msg = 'file={} incomplete header encountered'
            raise KeyError(msg.format(bfile))
        imax = int(self.info['CHROMOSOME_LENGTH'] / p_len)
        curr_pidx = self.info['CHR_PARTITION_INDEX']
        # update expected p_len for last segment in chrom
        if curr_pidx == imax:
            p_len = self.info['CHROMOSOME_LENGTH'] - (imax * p_len)
        # calculate the gap size
        gap_size = p_len - s_sum
        # fill gap in the last segment
        if gap_size <= 2:
            self._data[-1][1] += gap_size
            s_sum += gap_size
        else:
            msg = 'file={} segments do not sum to partition length'
            raise ValueError(msg.format(bfile))

        # safemode checks partition merge, has no effect for whole chrom maps
        if safemode:
            self._chkpart(last_pidx, s_sum)
        
    def write_bfile(self, f_out, safemode=True):
        """write a file header followed by the formatted bmap data"""
        # only write if b map data is complete
        if safemode:
            self._chkdata()
        # set pidx and plen to null values for merged map
        self._info['CHR_PARTITION_INDEX'] = -1
        self._info['CHR_PARTITION_LENGTH'] = -1
        with open(f_out, 'w') as f:
            for (key, value) in self.info.items():
                f.write('#{}={}\n'.format(key, value))
            for (b, s) in self.data:
                f.write('{:.0f} {:.0f}\n'.format(b, s))
            
    def merge_partitions(self, bfile_list):
        """process contiguous partitions and save to single file"""
        # re-initialize the reader and read the first file in list
        self.__init__(bfile=bfile_list[0])  # (this wipes data and info attrs)
        # read remaining files
        for bf in bfile_list[1:]:
            self.read_bfile(bf)
        
    def get_bmap(self, safemode=True):
        """returns a BkgdMap built from bfile data"""
        # only return if b map data is complete
        if safemode:
            self._chkdata()
        bmap = self.data_array
        segs = bmap[1]
        bkgd = self._b(bmap[0])
        return BkgdMap(bkgd, segs)
    
    def _hdrinfo(self, line):
        """format header info using regex, check contiguity of partitions"""
        key, value = line[1:].strip('\n').split('=')
        if re.match('^[+-]*\d+$', value):
            self._info[key] = int(value)
        elif re.match('^[-+]*\d+\.*[\d+-eE]+$', value):
            self._info[key] = float(value)
        else:
            self._info[key] = value

    def _b(self, b_int):
        """convert b integer bins to real-valued b values"""
        # get B and b precision parameters from info dict
        B_epsilon = 1.0 / self.info['BKGD_SCALE']
        b_epsilon = -np.log1p(-B_epsilon)
        return -b_int * b_epsilon
    
    def _chkpart(self, last_pidx, s_sum):
        """check for errors while merging partition maps"""
        # check partition contiguity (single partition files have pidx=-1)
        curr_pidx = self.info['CHR_PARTITION_INDEX']
        if (curr_pidx >= 0) and (curr_pidx != last_pidx + 1):
            msg = 'pidx={}: cannot join disjoint partition files'
            raise RuntimeError(msg.format(curr_pidx))
        # check that segments sum to expected partition length 
        p_len = self.info['CHR_PARTITION_LENGTH']
        imax = int(self.info['CHROMOSOME_LENGTH'] / p_len)
        if (imax > curr_pidx >= 0) and (s_sum != p_len):
            # msg = 'file={} segments do not sum to partition length'
            msg = 'pidx={}, len={}: segments do not sum to partition length'
            raise RuntimeError(msg.format(curr_pidx, s_sum))
        # check that the LAST segment sums to partition length
        last_plen = self.info['CHROMOSOME_LENGTH'] - (imax * p_len)
        if (curr_pidx == imax) and (s_sum < last_plen):
            msg = 'pidx={}: segments do not sum LAST partition length'
            raise RuntimeError(msg.format(curr_pidx))
    
    def _chkdata(self):
        """check if reader has loaded full chromosome-wide b map data"""
        if not self.data:
            raise RuntimeError('the reader has no data')
        if not self.info:
            raise RuntimeError('the reader is missing header info')
        if sum(d[1] for d in self.data) != self.info['CHROMOSOME_LENGTH']:
            msg = 'b map segments do not sum to {} length'
            raise RuntimeError(msg.format(self.info['CHROMOSOME_NAME']))
            
    @property
    def info(self):
        return self._info
    
    @property
    def data(self):
        return self._data
    
    @property
    def data_array(self):
        return np.array(self.data).T
