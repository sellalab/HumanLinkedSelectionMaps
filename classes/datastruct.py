from knowngene import KnownGene
from programs.calc_bkgd import calc_bkgd
from runstruct import chromosome_length, ChromStruct, izip, np, re
from wigfixblock import WigfixBlock, wig_iterator
from programs.compress_lsmaps import compress_lsmaps
from alignmentblock import AlignmentBlock, axt_iterator
from functions import count_all_pairs, relative_error, zopen

__author__ = 'davidmurphy'


def binary_mask_segments(mask):
    """
    find [start, end] indices for groups of 1 or more 1s in a binary {01}
    mask
    """
    # 1. create an array from the neutral mask with buffering 0s on either end
    zarr = np.concatenate(([0], mask, [0]))
    # 2. where mask != mask shifted by 1 delineate neutral segments edges
    edges = np.where(zarr[1:] != zarr[:-1])[0]
    # must be an equal number of starts and ends so len(edges) must be even
    assert len(edges) % 2 == 0
    # 3. reshape results (C-style) to get [start, end] coords of each segment
    edges = edges.reshape((len(edges) / 2), 2)

    return edges


def mask_segments(mask, segments, flipoff=True, returnmask=True):
    """
    Return mask with segments flipped to 0 or new segments after removing masked
    regions
    :param mask: chromosome spanning binary mask
    :param segments: 0-based annotation coordinates on the same chromosome
    :param flipoff: flip mask to False in segments if returning mask; map mask
    values onto segment regions if
    returning the segments (default). If flipoff is False, flip mask to True in
    segments if returning mask;
    map mask inverse onto segment regions when returning segments.
    :param returnmask: return mask after segment intersect when true (default),
    otherwise return new segments
    """
    # return the mask after altering regions in segments
    for (i, j) in segments:
        mask[i:j] = 0 if flipoff else 1
    return mask if returnmask else binary_mask_segments(mask)


def conservation_mask(fwig, chrom, pmin=0.0, pmax=0.0, mtype=bool, fout=None):
    """
    A function to read phastCons wiggle fixedstep formatted files and
    create a mask based on sites falling within a certain conservation
    score range. The mask is a simple binary mask with 1 for each position
    where the wigfix conservation file score is within [pmin, pmax]. The
    program can be used both with phastCons and phyloP files downloaded
    from UCSC. Missing sites are recorded 0s by default
    :param fwig: the wig file to use for getting cons/ncons regions
    :param chrom: name of masked chromosome
    :param pmin: miniumum acceptable score
    :param pmax: maximum acceptable score
    :param mtype: mask data type (e.g., bool, int, str)
    :param fout: optional file path to store finished mask
    """
    # set mask length
    mlen = chromosome_length(chrom)
    mask = np.zeros(mlen, dtype=np.bool)

    for bi in wig_iterator(wig_file=fwig):
        # make block objects from continuous score blocks in the wig file
        blk = WigfixBlock(data=bi.groups(), pmin=pmin, pmax=pmax)

        # make sure that blocks obey expected rules
        assert blk.chrom == chrom
        assert blk.step == 1

        # binarize scores into 0s and 1s according to pmin, pmax params
        mask[blk.start:blk.end] = blk.mask

    # update the mask to specified data type
    mask = mask.astype(mtype)

    # save to filename, if given
    if fout is not None:
        np.savez_compressed(fout, cmask=mask)

    return mask


def substitution_mask(faxt, chrom, fout=None):
    """
    Read an alignment file in .axt format, extract each block and then use
    the blocks to construct a mask file with the following score code:
        - 0: either no data on one or both sequences, or to pad a spatial
             gap between successive blocks
        - 1: true ACGT base in both sequences and seq1 == seq2
        - 2: true ACGT base in both sequences but seq1 != seq2

    :param faxt: optional faxt to process (default to stored
    filename in parent class)
    :param chrom: name of masked chromosome
    :param fout: a flag that indicates that the substitution mask
    should be saved to the default file location specified in RunStruct
    :return mask: numpy array the length of chrom using the 012 encoding
    described above or None if fout=True
    """
    # set mask length
    mlen = chromosome_length(chrom)
    mask = np.zeros(mlen, dtype='u1')

    # iterate through the alignment and build the mask
    num = -1  # block numbers given in the file
    end = 0  # current end pos
    for bi in axt_iterator(axt_file=faxt):
        # take the match object and make it into an AlignmentBlock object
        blk = AlignmentBlock(axt_block=bi.groups())

        # idx, pos increase correctly, seq sizes match, ch is correct
        assert blk.ref_chrom == chrom
        assert end <= blk.ref_start
        assert blk.alignment_number == num + 1
        assert len(blk.primary) == len(blk.aligned)

        # add block mask sequence at the block's chromosomal coordinates
        mask[blk.ref_start:blk.ref_end] = blk.submask()

        # reset positional indicators
        end = blk.ref_end
        num = blk.alignment_number

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, submask=mask)

    return mask


def neutrality_mask(nconsmask, callmask, submask, filter_files=(), fout=None):
    """
    A mask of the final set of fully filtered neutral positions
    :param nconsmask: mask of neutral, nonconserved sites -> can be created
    using conservation mask function
    :param callmask: mask of polymorphism calls
    :param submask: mask of match/mismatch with an outgroup
    :param filter_files: (optional) list of segment files to filter from mask
    :param fout: optional file path to store finished mask
    :return neutmask: FASTA-style mask encoded by 0s, 1s and 2s with the
    translated as follows:
        0=non-neutral/missing
        1=neutral & matches outgroup
        2=neutral & diverged from outgroup
    """
    # default mask: 100-vertebrates phastCons score == 0
    neutmask = nconsmask

    # remove segments from filter files
    for f in filter_files:
        segments = np.loadtxt(f, usecols=(1, 2), dtype='u4')
        mask_segments(neutmask, segments)

    # remove non-passing sites in polymorphism call mask (bool)
    neutmask &= callmask

    # remove sites without outgroup complement base for divergence
    neutmask = neutmask.astype('u1') * submask  # uint8 preserves 012 encoding

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, neutmask=neutmask)

    return neutmask


def neutral_snps(neutmask, snpcnt, fout=None):
    """
    Create an array of (0-based) [positions, sample_size, alt_count] for
    SNP data
    :param neutmask: FASTA-style mask encoded by 0s, 1s and 2s with the
    translated as follows:
        0=non-neutral/missing
        1=neutral & matches outgroup
        2=neutral & diverged from outgroup
    neutmask can be created with neutrality_mask function
    :param snpcnt: array of polymorphic sites formatted:
        SNP_POS REF_COUNT ALT_COUNT
    :param fout: optional file path to store finished mask
    :return snps:
    """
    # convert SNP positions to 0-based coordinate system
    pos = snpcnt[:, 0]
    pos -= 1  # NOTE: change numpy.view pos changes snpcnt as well

    # get boolean mask of neutral SNP sites
    isneut = neutmask[pos] > 0
    neutsnps = snpcnt[isneut]

    # save to filename, if given
    if fout:
        np.savez_compressed(fout, neutsnps=neutsnps)

    return neutmask


def snpcount(fsnp, returnbases=False):
    """
    Get polymorphism data from vcftools formatted .frq.count files
    :param fsnp: snp_count file formatted (after vcftools):
        CHROM POS N_ALLELES N_CHR REF:COUNT ALT:COUNT
    :param returnbases: optionally return the REF/ALT bases
    :return snp_array: [pos, refcount, altcount] or
    [pos, ref, refcount, alt, altcount] if returnbases=True
    ref, refcount, alt, altcount]
    """
    # match [snp_pos, ref_base, ref_count, alt_base, alt_count]
    if returnbases:
        re_snps = re.compile(
            '\d{1,2}\s(\d{1,9})\s\d+\s\d+\s([ACGT]):(\d+)\s([ACGT]):(\d+)'
        )
        with zopen(fsnp, 'r') as f:
            ar = np.array(re_snps.findall(f.read())).T
        snppos = ar[0].astype('u4')
        rbase = ar[1]
        rcount = ar[2].astype('u4')
        abase = ar[3]
        acount = ar[4].astype('u4')

        return snppos, rbase, rcount, abase, acount

    # match [snp_pos, ref_count, alt_count]
    else:
        re_snps = re.compile(
            '\d{1,2}\s(\d{1,9})\s\d+\s\d+\s[ACGT]:(\d+)\s[ACGT]:(\d+)'
        )
        with zopen(fsnp, 'r') as f:
            ar = np.array(re_snps.findall(f.read()), dtype='u4').T
        snppos = ar[0]
        rcount = ar[1]
        acount = ar[2]

        return snppos, rcount, acount


def sequence_mask(fcall, chrom, call_key='P'):
    """
    return a binary mask from calls, where passing=1, others=0
    :param fcall: path to sequencing call mask file
    :param chrom: name of masked chromosome
    :param call_key: symbol that indicates that a position passed filters
    (default, used by 1000 Genomes, is 'P')
    :return passed: a boolean array where called=True, error=False
    """
    # string -> list -> char array -> boolean mask: P=True, ~P=False
    with zopen(fcall, 'r') as m:
        _ = m.readline()  # skip header line
        chars = np.array(list(m.read().replace('\n', '')))
    assert chars.dtype == '|S1'
    assert chars.size == chromosome_length(chrom)
    passed = (chars == call_key)

    return passed


class DataStruct(ChromStruct):
    """Functions for extracting features from genome reference files"""
    def __init__(self, chrom='chr1', **kwargs):
        # TODO: any other intrinsic properties needed?
        super(DataStruct, self).__init__(chrom=chrom, **kwargs)
        # if 'clone' in kwargs:
        #     return
        self.chrom = chrom

    @property
    def refchrom(self):
        """return the reference for the current chromosome as a string array"""
        with zopen(self.refg_files, 'r') as f:
            f.readline()  # skip header
            ref = np.array(list(f.read().replace('\n', '').upper()))
        assert ref.dtype == '|S1'
        assert ref.size == self.chlen
        return ref

    @property
    def ancchrom(self):
        """return the reference for the current chromosome as a string array"""
        # NOTE: low conf calls are up-cased in this function!
        with zopen(self.ancs_files, 'r') as f:
            f.readline()  # skip header
            anc = np.array(list(f.read().replace('\n', '').upper()))
        assert anc.dtype == '|S1'
        assert anc.size == self.chlen
        return anc

    @property
    def conservation_percentiles(self):
        """return cutoff score for top 1-25% of cons score sites"""
        # look at the last 25-percentiles of the score distribution
        percentiles = np.arange(0.75, 1.0, 0.01)

        # convert counts to cumulative fraction of scores
        scores, counts = np.loadtxt(self.wigzcons_distribution).T
        counts = np.cumsum(counts) / counts.sum()

        # find the index of the first occurrence of each percentile
        idx = np.minimum(np.searchsorted(counts, percentiles), scores.size - 1)

        # return the corresponding p-scores mapped to their percentiles
        pct_keys = [str(int(p)) for p in 100*percentiles]
        pct2score = izip(pct_keys, scores[idx])
        return dict((p, s) for (p, s) in pct2score)

    @property
    def segments(self):
        """return segments from the prediction compression"""
        return np.concatenate([[0], np.cumsum(np.load(self.seg_files))])

    def callmask(self, mask_file=None):
        """
        return a binary mask from calls, where passing=1, others=0
        :param mask_file: specify path to non-default mask file
        """
        # string -> list -> char array -> boolean mask: P=True, ~P=False
        mask_file = self.call_mask if mask_file is None else mask_file
        with zopen(mask_file, 'r') as m:
            _ = m.readline()  # skip header line
            chars = np.array(list(m.read().replace('\n', '')))
        assert chars.dtype == '|S1'
        assert chars.size == self.chlen
        return chars == self.vars.callmask_key

    def substitution_mask(self, axt_file=None, save_file=False):
        """
        Read an alignment file in .axt format, extract each block and then use
        the blocks to construct a mask file with the following score code:
            - 0: either no data on one or both sequences, or to pad a spatial
                 gap between successive blocks
            - 1: true ACGT base in both sequences and seq1 == seq2
            - 2: true ACGT base in both sequences but seq1 != seq2
            
        :param axt_file: optional axt_file to process (default to stored
        filename in parent class)
        :param save_file: a flag that indicates that the substitution mask
        should be saved to the default file location specified in RunStruct
        :return mask: numpy array the length of chrom using the 012 encoding
        described above or None if save_file=True
        """
        # use internal values if no filepath is provided
        axt_file = self.axt_outgroup if axt_file is None else axt_file
        # empty zeros array to store the mask values returned from each block
        mask = np.zeros(shape=self.chlen, dtype='u1')

        # iterate through the alignment and build the mask
        num = -1  # block numbers given in the file
        end = 0  # current end pos
        for blk in axt_iterator(axt_file=axt_file):
            # take the match object and make it into an AlignmentBlock object
            b = AlignmentBlock(axt_block=blk.groups())

            # idx, pos increase correctly, seq sizes match, ch is correct
            assert b.ref_chrom == self.chrom
            assert end <= b.ref_start
            assert b.alignment_number == num + 1
            assert len(b.primary) == len(b.aligned)

            # add block mask sequence at the block's chromosomal coordinates
            mask[b.ref_start:b.ref_end] = b.submask()

            # reset positional indicators
            end = b.ref_end
            num = b.alignment_number

        if save_file:
            # save mask to axt_file path, replace extension w/ '.submask.npz'
            mask_file = self.axt_outgroup.replace('.net.axt.gz', '.submask.npz')
            np.savez_compressed(mask_file, submask=mask)
        else:
            return mask

    def substitution_counts(self, msk_file=None, window=None, save_file=False):
        """
        Record the number of match/mismatch neutral bases per specifed window
        size
        :param msk_file: saved mask file: open and count
        :param window: set different window size for counting neutral
        match/mismatch positions
        :param save_file: a flag that indicates that the substitution counts
        should be saved to the default file  location specified in RunStruct
        :return [window, bases, subs] array for the chrom
        If write flag set to True, files are written with the following format:
            # POS   BASES   SUBS
            100000  5432    834
        POS: upper bound of bin, i.e. rates from 0-1000 will have POS = 1000
        BASES: total number of "1" and "2" in the window, i.e. all valid bases
        SUBS: total number of "2" in the window, just the substitutions
        """
        # set window size to use for counting match/mismatch
        window = self.vars.muest_window if window is None else window
        # check if a filename is specified saved neutrality mask
        if msk_file is not None:
            neutmask = np.load(msk_file)['neutmask']
        # otherwise use default files to create neutrality mask
        else:
            neutmask = self.neutrality_mask()

        # calculate the total number of windows that will cover the chromosome
        tot_bins = int(np.ceil(1.0 * self.chlen / window))
        # empty array for [window_end, tot_bases, tot_subs] for each window
        counts_array = np.zeros(shape=(tot_bins, 3), dtype='u4')

        for idx in xrange(tot_bins):
            i, j = idx * window, (idx + 1) * window
            # count 0s, 1s and 2s in each window
            match, mismatch = np.bincount(neutmask[i:j], minlength=3)[1:]
            # record [upper_edge, tot_matches, tot_mismatches] per window
            counts_array[idx] = min(j, self.chlen), match + mismatch, mismatch

        # substitution rate estimates should tile the full chromosome
        assert counts_array[-1, 0] == self.chlen

        if save_file:
            np.savez_compressed(self.sub_count, subcount=counts_array)
        else:
            return counts_array

    def estimator_relative_error(self, msk_files=None, window=None):
        """
        Function to calculate the relative error in the mutation rate estimator
        across all autosomes
        :param msk_files: optional list of pre-processed neutrality mask files
        :param window: set different window size for counting neutral
        match/mismatch positions
        :return relative_err: a vector of the relative error in each window for
        all windows in the autosomes
        """
        msk_files = self.neut_masks if msk_files is None else msk_files
        window = self.vars.muest_window if window is None else window

        # length of each chrom in units of windows
        lengths = []
        for ch in self.chroms:
            nwin = int(np.ceil(1.0 * self.chromosome_length(ch) / window))
            lengths.append(nwin)

        # record [bases, subs] for each bin in every chrom
        counts_array = np.zeros(shape=(sum(lengths), 2))
        i = 0
        for idx in xrange(22):
            self.chrom = self.chroms[idx]
            # results/chrom[i] fill exactly lengths[i] rows in counting array
            counts = self.substitution_counts(msk_files[idx], window)[:, 1:]
            counts_array[i:i + lengths[idx]] = counts
            i += lengths[idx]

        # mask low data windows, calculate relative error for remaining data
        bases, subs = counts_array.T
        mask = (bases >= self.vars.min_bases)
        bases, subs = bases[mask], subs[mask]
        err_array = relative_error(n=bases, k=subs)

        # print the mean relative error
        kb = '{}kb'.format(int(window / 1000))
        vals = kb, self.vars.min_bases, err_array.mean(), np.median(err_array)
        print 'window={}; minbases={}; mean={}; median={}'.format(*vals)

        return err_array

    # todo: define upper/lower scores, i.e. neutral~[0, x] vs neutral~[0,0]
    def conservation_mask(self, fwig, pct=None, masktype=bool, save_file=False):
        """
        A function to read phastCons wiggle fixedstep formatted files and
        create a mask based on sites falling within a certain conservation
        score range. The mask is a simple binary mask with 1 for each position
        where the wigfix conservation file score is within [pmin, pmax]. The
        program can be used both with phastCons and phyloP files downloaded
        from UCSC. Missing sites are recorded 0s by default
        :param fwig: the wig file to use for getting cons/ncons regions
        :param masktype: the data type used to record from the mask
        (default is bool)
        :param pct: the percent cutoff to use for calling a site
        conserved (e.g., 0.95 will return only sites with conservation scores
        in the top 95th percentile). use 0 for neutral only (scores=0)
        :param save_file: a flag to save conserved segments to the bed_file
        specified by RunStruct
        """
        # fully neutral case
        if pct == 0:
            pmin = pmax = 0.0
        # custom percentile specified
        elif pct > 0:
            # convert percentile conservation percentile key (.2f string)
            pmin, pmax = self.conservation_percentiles[str(pct)], 1.0
        # default percentile used
        else:
            default_pct = self.vars.percentile
            pmin, pmax = self.conservation_percentiles[default_pct], 1.0

        # initialize appropriate mask type
        if masktype is bool:
            mask = np.zeros(shape=self.chlen, dtype=np.bool)
        elif masktype is np.array:
            mask = np.full(shape=self.chlen, fill_value=np.nan, dtype='u2')
        elif masktype is np.uint8:
            mask = np.zeros(shape=self.chlen, dtype='u1')
        else:
            mask = np.zeros(shape=self.chlen, dtype='S1')

        # generate blocks from wig file
        for blk in wig_iterator(wig_file=fwig):
            # expand the match object into a WigFixBlock object
            b = WigfixBlock(data=blk.groups(), pmin=pmin, pmax=pmax)

            # make sure that the block obeys the normal rules
            assert b.chrom == self.chrom
            assert b.step == 1

            # record scores in appropriate data type for the current block
            if masktype is bool:
                mask[b.start:b.end] = b.mask_bool
            elif masktype is np.array:
                mask[b.start:b.end] = b.score_array
            elif masktype is np.uint8:
                mask[b.start:b.end] = b.mask_int
            else:
                mask[b.start:b.end] = b.mask_str

        if save_file:
            # save bed to fwig path, change extension to '.pct.bed'
            suff = '{}.bed'.format(pct)
            bed_file = fwig.replace('.wig.gz', suff)
            # convert mask to segments
            segments = binary_mask_segments(mask)
            # create column of "chrN" for bed formatted file
            stype = 'S{}'.format(len(self.chrom))
            column_1 = np.full(len(segments), self.chrom, dtype=stype)
            # join chrom label, start, end coords into bed array
            bed_array = np.column_stack((column_1, segments.astype(str)))
            np.savetxt(bed_file, bed_array, fmt='%s\t%s\t%s')
        else:
            return mask

    def neutrality_mask(self, wigmask=None, callmask=None, submask=None,
                        save_file=False):
        """
        A mask of the final set of fully filtered neutral positions
        :param wigmask: mask of neutral, nonconserved sites
        :param callmask: mask of polymorphism calls
        :param submask: mask of match/mismatch with an outgroup
        :param save_file: a flag to save the mask file to default file location
        :return neutrality_mask: 012-encoded mask of neutral sites,
        where 0=masked, 1=neutral match, 2=neutral mismatch (vs. outgroup)
        """
        # count the remaining neutral sites at each filtering step
        site_counter = [('initial', self.chlen)]

        # default mask: 100-vertebrates phastCons score == 0
        if wigmask is None:
            mask = self.conservation_mask(self.wigz_ncons, pct=0.0)
        # use wigmask from args
        else:
            mask = wigmask
        site_counter.append(('phastCons', np.sum(mask != 0)))

        # remove sites overlapping bs targets (i.e., conserved regions)
        if self.vars.filter_con:
            for f in self.bs_targets.values():
                mask_segments(mask, np.loadtxt(f, usecols=(1, 2), dtype='u4'))
        site_counter.append(('conserved', np.sum(mask != 0)))

        # remove positions in filtered regions (e.g., genes, duplications...)
        for f in self.filter_files.values():
            mask_segments(mask, np.loadtxt(f, usecols=(1, 2), dtype='u4'))
        site_counter.append(('filtered', np.sum(mask != 0)))

        # remove non-passing sites in polymorphism call mask (bool)
        mask &= self.callmask() if callmask is None else callmask
        site_counter.append(('callmask', np.sum(mask != 0)))

        # remove sites without outgroup complement base for divergence
        submask = self.substitution_mask() if submask is None else submask
        # use uint8 for (mask) x (submask) to preserve 012 divergence encoding
        mask = mask.astype('u1') * submask
        site_counter.append(('submask', np.sum(mask != 0)))

        if save_file:
            np.savez_compressed(self.neut_masks, neutmask=mask)
        else:
            return mask

    def neutral_snps(self, msk_file=None, snpcount=None, save_file=False):
        """
        Create an array of (0-based) [positions, sample_size, alt_count] for
        SNP data
        :param msk_file: FASTA style mask coded 012 with 1s are neutral sites
        matching the outgroup and 2s are neutral sites diverged from outgroup.
        In the final output divergent sites are represented by 1s and matching
        sites by 0s (default neutrality_mask generated from scratch using files
        specified in the parent class)
        :param snpcount: array of polymorphic sites including refcount and
        altcount data (default array generated from stored files)
        :param save_file: flag to specify that array should be saved to
        default file location
        """
        # prepare neutrality_mask
        msk_file = self.neut_masks if msk_file is None else msk_file
        neutmask = np.load(msk_file)['neutmask']

        # get SNP data, use 0-based to match the neutrality_mask coordinates
        spos, ref, alt = self.snpcount() if snpcount is None else snpcount
        spos -= 1

        # keep SNPs at neutral sites
        isneut = neutmask[spos] > 0
        snps = np.column_stack((spos[isneut], ref[isneut], alt[isneut]))

        if save_file:
            np.savez_compressed(self.neut_poly, neutpoly=snps)
        else:
            return snps

    def snpcount(self, fsnp=None, returnbases=False):
        """
        Get polymorphism data from vcftools formatted .frq.count files
        :param fsnp: snp_count file formatted (after vcftools):
            CHROM POS N_ALLELES N_CHR REF:COUNT ALT:COUNT
        :param returnbases: optionally return the REF/ALT bases
        :return snp_array: [pos, refcount, altcount] or
        [pos, ref, refcount, alt, altcount] if returnbases=True
        """
        fsnp = self.snp_files if fsnp is None else fsnp

        # match [snp_pos, ref_base, ref_count, alt_base, alt_count]
        if returnbases:
            re_snps = re.compile('\d{1,2}\s(\d{1,9})\s\d+\s\d+\s'
                                 '([ACGT]):(\d+)\s([ACGT]):(\d+)')
            with zopen(fsnp, 'r') as f:
                ar = np.array(re_snps.findall(f.read())).T
                snppos = ar[0].astype('u4')
                rbase = ar[1]
                rcount = ar[2].astype('u4')
                abase = ar[3]
                acount = ar[4].astype('u4')
            return snppos, rbase, rcount, abase, acount

        # match [snp_pos, ref_count, alt_count]
        else:
            re_snps = re.compile('\d{1,2}\s(\d{1,9})\s\d+\s\d+\s'
                                 '[ACGT]:(\d+)\s[ACGT]:(\d+)')
            with zopen(fsnp, 'r') as f:
                ar = np.array(re_snps.findall(f.read()), dtype='u4').T
                snppos = ar[0]
                rcount = ar[1]
                acount = ar[2]
            return snppos, rcount, acount

    def cpgsnps(self, fsnp=None, save_file=False):
        """return the positions of SNPs resulting from CpG
        transitions relative to the reference"""
        # get neutral SNP positions and base info
        snppos, ref_base, _, alt_base, _ = self.snpcount(fsnp, True)

        # get positions of G>A transitions and C>T pos
        gapos = snppos[((ref_base == 'G') & (alt_base == 'A'))]
        ctpos = snppos[((ref_base == 'C') & (alt_base == 'T'))]

        # keep the positions of confirmed CpG transitions: CG>CA and CG>TG
        ref_chrom = self.refchrom
        # -1 for 0-based index and -1 to get 5' flanking site of G>A
        ga_cpg = gapos[ref_chrom[gapos - 2] == 'C']
        # -1 for 0-based index and +1 to get 3' flanking site of C>T
        ct_cpg = ctpos[ref_chrom[ctpos] == 'G']

        # concatenate and return or save the sorted CpG positions
        if save_file:
            pass
        else:
            return np.sort(np.concatenate((ga_cpg, ct_cpg)))

    def conserved_annotations(self, conserved=None, annotations=('exons',)):
        """
        Return each of the annotations specified by keyword args partitioned by
        overlap with conserved regions
        :param conserved: optional path to conserved segments (use RunStruct
        specified regions by default)
        :param annotations: list of elements from KnownGene class
        """
        # read input conservation scores file or default to store file
        conserved_file = self.bed_cons if conserved is None else conserved

        # create a mask of conserved regions
        cmask = np.zeros(shape=self.chlen, dtype='u1')
        for (i, j) in np.loadtxt(conserved_file, usecols=(1, 2), dtype='u4'):
            cmask[i:j] = 1

        # load genes as KnownGenes for the chrom
        with open(self.genes, 'r') as f:
            for line in f:
                g = KnownGene(line)
                for anno in annotations:
                    # multiply annotation segments by 2
                    for (i, j) in getattr(g, anno):
                        # (will only affect conserved since neutral=0)
                        cmask[i:j] *= 2

        # conserved in annotation
        cons_in = binary_mask_segments(mask=(cmask > 1))

        # conserved outside annotation
        cons_out = binary_mask_segments(mask=(cmask == 1))

        return cons_in, cons_out

    """
    PRE-CALCULATION AND COMPRESSION FUNCTIONS
    """
    def calc_bkgd(self, tval, bsanno, save_cfg=False):
        """
        Call McVicker's calc_bkgd program on each bs_anno for the deleterious
        fitness effect "tval"
        :param tval: the deleterious fitness effect (mean should be >> 2Ns)
        :param bsanno: background selection target annotation
        :param save_cfg: flag to save the config file used with calc_bkgd
        """
        if isinstance(tval, str):
            tval = eval(tval)
        else:
            assert isinstance(tval, float)
        calc_bkgd(cst=self, tval=tval, bsanno=bsanno, save_cfg=save_cfg)

    def compress_predictions(self):
        """
        The purpose of this script is to 'digitize' McVicker bmaps so that
        for each position where one of a set of maps changes to a new B
        value, the other maps are 'cut' at this site as well. The end
        product should be a set of maps for different selection coefficients
        but with the exact same coordinates for each 'cut'. These digitized
        maps can then be used to squeeze the maximum number of neutral sites
        into each segment from the set of bmaps. We can then pool the
        neutral data that fit into a given bmap segment to avoid calculating
        the likelihood at every single site at each iteration of the
        composite likelihood calculation. We also build the mutation rate
        estimate into the minimum set of segments where estimates are
        available.
        """
        bmaps = []
        cmaps = []

        # flatten bs/cs into single lists
        for anno in self.bs_annos:
            bmaps.extend(self.bmap_files[anno])
        for anno in self.cs_annos:
            cmaps.extend(self.cmap_files[anno])
        subs = self.substitution_counts(msk_file=self.neut_masks)

        # run compression program
        seg, bs, cs, mu = compress_lsmaps(bmaps, cmaps, subs, self.chlen)

        # save outputs
        np.save(self.seg_files, seg)
        np.save(self.u_files, mu)
        if bmaps:
            np.save(self.bs_files, bs)
        if cmaps:
            np.save(self.cs_files, cs)

    def compress_data(self, msk_file=None):
        """
        compresses polymorphism data array by grouping sites into discrete
        prediction windows
        """
        # get neutral data arrays ready
        msk_file = self.neut_masks if msk_file is None else msk_file
        neutmask = np.load(msk_file)['neutmask']
        snps = np.load(self.neut_poly)['neutpoly']

        # check that every position in SNPs has the same sample size
        sample = np.sum(snps[:, 1:], axis=1)
        assert sample.min() == sample.max()

        # use this for calculations below
        sample = sample[0]
        combos = 0.5 * (sample**2 - sample)

        # use histogram on SNP positions to get SNP site indices per seg idx
        numbins = np.maximum(0, self.segments-1)
        sidx, nidx = np.histogram(a=snps[:, 0], bins=numbins)
        sidx = np.cumsum(np.concatenate(([0], sidx)))
        n = len(sidx) - 1  # number of segments

        # array for data summaries in each segment
        summaries = np.zeros(shape=(n, 4), dtype='f8')

        for i in xrange(n):
            # indices of the SNP data based on counts
            si, sj = sidx[i:i+2]
            polymorph = sj - si
            # sum the homo/het pairs if segment contains any SNPs
            if polymorph:
                alt = snps[si:sj, -1].T
                hom, het = map(np.sum, count_all_pairs(sample, alt))
            else:
                hom = het = 0.0

            # neutral mask indices are just segment boundaries
            ni, nj = nidx[i:i+2]
            # use bincount to sum match=1, mismatch=2 sites in the segment
            match, mismatch = np.bincount(neutmask[ni:nj], minlength=3)[1:]
            total = match + mismatch

            # count homozygous combos for non-SNP neutral sites in the segment
            if total:
                hom += (combos * (total - polymorph))

            # record summary information
            pairs = hom + het
            assert (pairs / combos == total) and (pairs % combos == 0)
            summaries[i] = hom, het, total, mismatch

        np.save(self.nt_files, summaries[:, :2])
        np.save(self.ndiv_files, summaries[:, 2:])
