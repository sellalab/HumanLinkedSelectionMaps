from mapstruct import RunStruct, MapStruct, np
from functions import chromosome_length, count_all_pairs, os

__author__ = 'davidmurphy'


class HiPiData(MapStruct):
    """
    A subclass of the MapStruct for analyzing high diversity data
    """
    def __init__(self, frac=0.02, **kwargs):
        """
        Saves the index in masked data for the top 'pct' % of prediction sorted data segments.
        :param frac: the fraction of highest diversity prediction data to take
        """
        super(HiPiData, self).__init__(**kwargs)
        pct = 100 * frac
        # data directory
        self.data_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/data'
        # the output file name for high diversity data segment indices
        self.hipi_file = '{}/hipi/{}.{}.{}pct.hipi.idx.npy'.format(self.root, self.neut, self.token, int(pct))
        # high diversity polymorphic sites file
        self.hipi_poly = '{}/hipi/{}.{}pct.hipi.snp.npy'.format(self.root, self.token, int(pct))
        # high diversity cpg indices file
        self.hipi_cpgs = '{}/hipi/{}.{}pct.hipi.rcpg.idx.npy'.format(self.root, self.token, int(pct))
        # all data cpg indices file
        self.all_cpgs = '{}/snps/{}.rcpg.idx.npy'.format(self.root, self.token)

        # open if data already saved, otherwise generate from MapStruct
        if os.path.isfile(self.hipi_file):
            self.hidx = np.load(self.hipi_file)
        else:
            self.hidx = self.hipi_index(p=frac)
            np.save(file=self.hipi_file, arr=self.hidx)

    def hipi_index(self, p=0.02):
        """
        Use prediction sort indices to select data with the highest predicted diversity. Keep the indices of the upper 
        'p' fraction of these high diversity data.
        :param p: fraction of top sites to keep - default = 2%
        :return hipi_index: indices corresponding to the top 'p' high diversity sites in masked dataset
        """
        sort_idx = self.prediction_sorting_index  # slice past sort_idx - no longer need it
        data_idx = np.arange(self.mcnt)[sort_idx]  # sort data indices by predicted diversity
        cum_sites = np.cumsum(self.nsites[sort_idx])  # take cumulative sum of sorted counts and keep final 'p' sites
        top_idx = np.where(cum_sites > (1 - p) * cum_sites.max())[0]
        return data_idx[top_idx]

    def compare_feature(self, feature, remove_self=False, random_sample=False):
        """
        Compare high diversity sites vs. genome (or genome - high diversity) sites by the fraction of neutral sites
        falling in segments with high vs. low fraction of some feature (e.g., callability mask). Generate a histogram
        based on the counts of neutral sites in a range from 0-100% passing for the two data sets.
        :param feature: a feature of the segments to compare for high diversity and genomic segments
        :param remove_self: optional flag to remove the high diversity segments from 'all' segments
        :param random_sample: replace high diversity with a random sample of the name number of segments from anywhere
        :return hdhist, gnhist: histogram counts for high diversity and genome sites counting sites in % feature bins
        """

        # get the masked feature, segments and nsites from genomic data
        gnfeat, gnseg, gnsites = self.masked([feature, 'seg', 'nsites'])
        # feature fraction in each genomic segment
        gnfeature = gnfeat.astype(float) / gnseg.astype(float)

        if random_sample:
            # take a random set of data segments to compare with the total set
            ridx = np.random.randint(low=0, high=len(gnfeature), size=len(self.nsites))
            hdfeature = gnfeature[ridx]
            hdnsites = gnsites[ridx]
        else:
            # take the high diversity sites using the high diversity index
            hdfeature = gnfeature[self.hidx]
            hdnsites = gnsites[self.hidx]

        if remove_self:
            # remove high diversity or random_sample feature fraction and site counts from total sites
            gnfeature = np.delete(gnfeature, self.hidx)
            gnsites = np.delete(gnsites, self.hidx)

        # get histogram counts for each set of features
        hdhist = np.histogram(hdfeature, bins=100, weights=hdnsites)
        gnhist = np.histogram(gnfeature, bins=100, weights=gnsites)

        return hdhist[0], gnhist[0]

    def compare_repeat(self, mst, repeat, remove_self=False, random_sample=False):
        """
        Compare the fraction neutral sites in high diversity and genome segments with different fractions of repeats
        present.
        :param repeat: the repeat class to look at
        :param remove_self: optional flag to remove the high diversity segments from 'all' segments
        :param random_sample: replace high diversity with a random sample of the name number of segments from anywhere
        :type mst: MapStruct
        :return hdhist, gnhist: histogram counts for high diversity and genome sites counting sites in % feature bins
        """

        # empty lists to accumulate repeats into
        hdr = []
        gnr = []
        gnranges = mst.masked(['ranges'])[0]
        segments = mst.masked(['seg'])[0]

        if random_sample:
            ridx = np.random.randint(low=0, high=len(gnranges), size=len(self.nsites))
            hdranges = gnranges[ridx]
            hdchroms = mst.masked(['chrom_tag'])[0][ridx]
            hdnsites = mst.masked(['nsites'])[0][ridx]
        else:
            hdranges = gnranges[self.hidx]
            hdchroms = mst.masked(['chrom_tag'])[0][self.hidx]
            hdnsites = mst.masked(['nsites'])[0][self.hidx]

        # load repeat array file if it already exists
        # repeat_array = '{dr}/{rp}/{lb}_repeat_array.npy'.format(dr=self.rep_dir, rp=repeat, lb=mst.label)
        repeat_array = '{dr}/coords/{rp}/{rp}_array.npy'.format(dr=self.data_dir, rp=repeat)
        if os.path.isfile(repeat_array):
            rpa = np.load(repeat_array)[mst.ms]
        # otherwise create from each chrom
        else:
            rpa = self._reparray(mst=mst, repeat=repeat)[mst.ms]

        # count repeats in segments for high diversity and genome data for each chrom
        for ch in mst.chroms:
            # get the segments ranges for the current chrom for high diversity and genome data
            cidx = mst.chroms.index(ch)
            msi, msj = mst.msij[cidx]
            hdrange = hdranges[hdchroms == cidx + 1]
            gnrange = gnranges[msi:msj]
            segs = segments[msi:msj]
            # subset of site indices from high diversity
            hidx = np.in1d(gnrange[:, 0], hdrange[:, 0])
            chrom_array = rpa[msi:msj]
            # now get fraction of repeats in high diversity and genome segments
            hdr.append([float(x) / float(s) for (x, s) in zip(chrom_array[hidx], segs[hidx])])
            gnr.append([float(x) / float(s) for (x, s) in zip(chrom_array, segs)])

        # get histogram counts for each set of repeats
        hdhist = np.histogram(np.concatenate(hdr), bins=100, weights=hdnsites)
        gnhist = np.histogram(np.concatenate(gnr), bins=100, weights=mst.masked(['nsites'])[0])

        return hdhist[0], gnhist[0]

    def _reparray(self, mst, repeat):
        """generate and save repeat array for repeat type specified"""
        # repeat_array = '{dr}/{rp}/{lb}_repeat_array.npy'.format(dr=self.rep_dir, rp=repeat, lb=mst.label)
        repeat_array = '{dr}/coords/{rp}/{rp}_array.npy'.format(dr=self.data_dir, rp=repeat)
        rpa = []
        chrom_ranges = mst.ranges
        for ch in mst.chroms:
            # repeat_file = '{dr}/{rp}/{ch}.{rp}.hg19.bed.gz'.format(dr=self.rep_dir, rp=repeat, ch=ch)
            repeat_file = '{dr}/coords/{rp}/{ch}.{rp}.bed'.format(dr=self.data_dir, rp=repeat, ch=ch)
            # get the segment ranges for the current file
            cidx = mst.chroms.index(ch)
            ci, cj = mst.ij[cidx]

            # start with an empty zeros the length of the chrom, flip to 1s at repeat segments
            chrom_array = np.zeros(shape=chromosome_length(ch), dtype='u1')
            # with zopen(repeat_file, 'r')as f:
            with open(repeat_file, 'r')as f:
                f.next()  # skip header
                for line in f:
                    i, j = map(int, line.split('\t')[1:3])
                    chrom_array[i:j] = 1
            # count the repeat bases at each segment
            rpa.append([np.sum(chrom_array[i:j], dtype='f8') for (i, j) in chrom_ranges[ci:cj]])

        # concatenate counts across chroms
        rpa = np.concatenate(rpa)
        # save to file and return the array
        np.save(file=repeat_array, arr=rpa)

        return rpa

    @property
    def snps(self):
        """return [ch, pos, sample, alt] for all high diversity SNP sites"""
        if os.path.isfile(self.hipi_poly):
            return np.load(self.hipi_poly)
        else:
            return self._snps()

    def _snps(self):
        """generate the [ch, pos, sample, alt] array for all high diversity SNP sites, save and then return"""
        # key rows for decompressing data segments to find original SNP sites
        hi = self.hidx
        # create array of cumulative sites per chrom
        cumsites = np.concatenate([np.cumsum(self.nsites[i:j]) for (i, j) in zip(self.mbnds[:-1], self.mbnds[1:])])
        # join arrays together into one array for common masking/sorting operations
        vs = np.column_stack([ar[hi] for ar in self.chrom_tag, cumsites, self.nsites, self.nt])
        hdd_poly = []
        for c in range(1, 23):
            # isolate rows for current chrom
            ar = vs[vs[:, 0] == c]
            # sort rows by cumulative neutral site count
            ar = ar[ar[:, 1].argsort()]
            # neutral data array start positions = cumsites - nsites for current seg
            si = ar[:, 1] - ar[:, 2]
            # neutral data array end positions = cumsites
            sj = ar[:, 1]

            # expand the individual sites from within each high diversity data segment using start, end coords
            row = 0
            chrom_poly = []
            neut = np.load(self.neut_poly[c - 1])['neutpoly']
            for (i, j) in zip(si, sj):
                # save all polymorphic sites to current list
                sample, alt = neut[i:j, 1:3].T
                chrom_poly.append(neut[i:j][np.where(alt > 0)[0]])

                # sanity check 1: reconstructed segments have the same het, hom counts
                hom, het = count_all_pairs(sample, alt)
                assert np.sum(hom) == ar[row, 3]
                assert np.sum(het) == ar[row, 4]

                row += 1

            # concatenate poly site data for current chrom
            hdd_poly.append(np.concatenate(chrom_poly))

        # create an array of concatenated data with one entry per chrom and save, then return
        final_array = np.array(hdd_poly)
        np.save(file=self.hipi_poly, arr=final_array)

        return final_array

    @property
    def cpgs_ref(self):
        """return [ch, pos, sample, alt] for all high diversity CpG SNPs polarized by reference"""
        if os.path.isfile(self.hipi_cpgs):
            idx = np.load(self.hipi_cpgs)
            return [a[i] for (a, i) in zip(self.snps, idx)]
        else:
            idx = self._cpgidx(subset=self.snps, subfile=self.hipi_cpgs)
            return [a[i] for (a, i) in zip(self.snps, idx)]

    @property
    def cpgs_anc(self):
        """return [ch, pos, sample, alt] for all high diversity CpG SNPs polarized by ancestor"""
        acpg_file = self.hipi_cpgs.replace('rcpg_idx', 'acpg_idx')
        if os.path.isfile(acpg_file):
            idx = np.load(acpg_file)
            return [a[i] for (a, i) in zip(self.snps, idx)]
        else:
            idx = self._cpgidx(subset=self.snps, subfile=acpg_file, ancestor=True)
            return [a[i] for (a, i) in zip(self.snps, idx)]
