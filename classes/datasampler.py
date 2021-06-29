from runstruct import RunStruct, np, os
from functions import chromosome_length

__author__ = 'davidmurphy'


class DataSampler(RunStruct):
    """
    A simple class to manage bootstrapping/LOOCV/Jackknife sampling of genomic data "in place" via array masking. The
    DataSampler class subclasses RunStruct, since it uses many of the same properties
    """
    def __init__(self, **kwargs):
        super(DataSampler, self).__init__(**kwargs)
        self.sample_type = self.vars.sample_type
        self.idx = []
        self.bootstrap_vector = []
        self.jackknife_vector = []

    def get_mask(self):
        """return the sampling mask vector"""
        # if no sampling type is specified by the RunStruct then no mask is returned
        if self.sample_type is None:
            return None
        if self.sample_type == 'boot':
            sample = self.get_boot()
        elif self.sample_type == 'jack':
            sample = self.get_jack()
        else:
            raise ValueError('unrecognized sample type "{}"'.format(self.sample_type))
        indices = self.get_idx()
        # start with all 0 weights for every position in the data
        coefficient_mask = np.zeros(shape=indices[-1, 1])
        # for each block region sampled, increase the mask by 1
        for (start, end) in indices[sample, :2]:
            coefficient_mask[start:end] += 1
        return coefficient_mask

    def get_idx(self):
        """return indices partitioning data into blocks"""
        # if another process has already generated the index no need to make it again
        if len(self.idx):
            return self.idx
        # check if the indices file exists
        elif os.path.isfile(self.data_indices):
            self.idx = np.load(self.data_indices)
        # otherwise create and save
        else:
            self.idx = self.idx_array()
            np.save(file=self.data_indices, arr=self.idx)
            # record the number of of total indices in the jackknife if it is not yet set
            if self.vars.jackknife_samples != np.sum(self.idx[:, 2]):
                self.vars.jackknife_samples = np.sum(self.idx[:, 2])
                self.save()

        return self.idx

    def get_boot(self, new=False):
        """
        return stored bootstrap sample or generate new
        :param new: flag indicates that a new sample should be created at the current index
        """
        # if internal bootstrap sample exists just use this
        if len(self.bootstrap_vector):
            return self.bootstrap_vector
        # check for existing bootstrap sample and check that new_sample is set to false, otherwise make new and save
        elif os.path.isfile(self.bootstrap_file) and not new:
            self.bootstrap_vector = np.load(self.bootstrap_file)
        else:
            idx = self.get_idx()
            # make a new temp dir for new bootstrap samples (incorporate the yy-mm-dd date in the name)
            newdir = '{}/boot/new_samples_{}'.format(self.root, self.record_timestamp()[:6])
            if not os.path.isdir(newdir):
                os.mkdir(newdir)
            newfile = '{}/sample-{}.npy'.format(newdir, self.vars.bootstrap_iteration)
            # use the final column to mask rows rows with no data
            hasdata = idx[:, 2].astype(bool)
            # the sample space only includes blocks which contain data
            sample_space = np.arange(0, len(idx))[hasdata]
            # the number of samples equals the number of blocks containing data in the original sample
            self.bootstrap_vector = np.random.choice(a=sample_space, size=len(sample_space), replace=True)
            # save sampled indices for re-use
            np.save(newfile, arr=self.bootstrap_vector)

        return self.bootstrap_vector

    def get_jack(self):
        """return stored jackknife sample or generate new"""
        # recycle internal variable
        if len(self.jackknife_vector):
            return self.jackknife_vector
        else:
            # create a jackknife sample of the data indices
            idx = self.get_idx()
            # use the final column to mask rows rows with no data
            hasdata = idx[:, 2].astype(bool)
            # the sample space only includes blocks which contain data
            sample_space = np.arange(0, len(idx))[hasdata]
            # pop out the jackknife block and save the sample to internal variable
            self.jackknife_vector = np.delete(sample_space, self.vars.jackknife_index)
        return self.jackknife_vector

    def idx_array(self):
        """create an array of [start, end, flag] rows for block ranges with flag=0 for empty blocks"""
        idx, cnts = self.split_data()
        start = 0
        idx_data = []
        for end in idx:
            if np.sum(cnts[start:end]):
                # indices for sample regions containing ANY neutral data (final column is a boolean flag)
                idx_data.append([start, end, 1])
            else:
                # indices for sample regions with NO data
                idx_data.append([start, end, 0])
            start = end
        return np.array(idx_data)

    def split_data(self):
        """get the segments (should tile the full length of the 22 autosomes)"""
        segs = np.concatenate([np.load(f) for f in self.seg_files])
        # get the cumulative positions across autosomes
        cumpos = np.cumsum(segs)
        # check that the segments tile the complete 22 autosomes (or chroms in the RunStruct)
        assert cumpos[-1] == sum(chromosome_length(chrom=ch) for ch in self.chroms)
        # get the site counts used for divergence calculations for each segment (0s are segments containing no data)
        cnts = np.concatenate([np.load(f)[:, 0] for f in self.ndiv_files])
        assert len(segs) == len(cnts)
        # find the indices delimiting non-overlapping 'sample_region' groups of segments
        sample_region = self.vars.sample_window
        idx = np.searchsorted(a=cumpos, v=np.arange(sample_region, cumpos[-1] + sample_region, sample_region))
        # return the indices separating sample blocks and the counts of neutral positions in each block
        return idx, cnts

    def next_jack(self):
        """move to the next jackknife index and reset the internal vector variable"""
        self.vars.jackknife_index += 1
        self.jackknife_vector = []

    def next_boot(self):
        """move to the next bootstrap index and reset the internal vector variable"""
        self.vars.bootstrap_iteration += 1
        self.bootstrap_vector = []

    @property
    def bootstrap_file(self):
        """get name of bootstrap indices file for current bootstrap iteration (it may not actually exist yet)"""
        sample_folder = '{}/boot/samples'.format(self.root)
        return '{}/sample-{}.npy'.format(sample_folder, self.vars.bootstrap_index)
