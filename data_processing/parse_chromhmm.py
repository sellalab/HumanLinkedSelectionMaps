__author__ = 'davidmurphy'


import os
import numpy as np
from collections import defaultdict
from data_processing.merge import merge
from data_processing.data_tools import binary_mask_segments
from classes.runstruct import root_dir, human_autosomes, chromosome_length
coords_dir = root_dir + '/data/coords'


class ChromHMMEntry:
    def __init__(self, data_string):
        """get elements of chromHMM from a string of data"""
        data = data_string.split()
        self.chrom = data[0]
        self.start = int(data[1])
        self.end = int(data[2])
        self.name = data[3].strip()


class ChromHMMDatabase:
    """a data structure for genome wide chromHMM ccordinates and labels"""
    anno_dict = dict(txn=('9_Txn_Transition', '10_Txn_Elongation'),
                     pro=('1_Active_Promoter', '2_Weak_Promoter',
                          '3_Poised_Promoter'),
                     enh=('4_Strong_Enhancer', '5_Strong_Enhancer'),
                     ins=('8_Insulator',),
                     rep=('12_Repressed',))
                          # '6_Weak_Enhancer', '7_Weak_Enhancer','8_Insulator'))

    def __init__(self, cell_line):
        """load and sort data by chroms and chromHMM names"""
        # format path to data files
        f_fmt = coords_dir + '/wgEncodeBroadHmm{}HMM.bed'
        data_file = f_fmt.format(cell_line)
        self.cell_line = cell_line
        # create nested dict for data by chrom and then by name
        self.chromdict = defaultdict(lambda: defaultdict(list))
        with open(data_file, 'r') as f:
            for line in f:
                hmm = ChromHMMEntry(line)
                self.chromdict[hmm.chrom][hmm.name].append(hmm)

    def _mkdir(self, anno):
        """make and return save path"""
        d_fmt = root_dir + '/data/bsanno/{}_{}'
        out_dir = d_fmt.format(self.cell_line.lower(), anno)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        return out_dir

    def _savefile(self, anno, out_dir, chrom):
        """get path to save file"""
        f_fmt = '{}/{}.{}_{}.bed'
        return f_fmt.format(out_dir, chrom, self.cell_line.lower(), anno)

    def get_chrom(self, chrom, name):
        """get single data type for a chromosome"""
        return self.chromdict[chrom][name]

    def get_segments(self, chrom, name):
        """get segments of a certain type for a given chrom"""
        return [(hmm.start, hmm.end) for hmm in self.chromdict[chrom][name]]

    def create_anno(self, anno):
        """save select annotation for all chroms"""
        tot = 0
        out_dir = self._mkdir(anno)
        for ch in human_autosomes:
            segs = []
            # get all the segments for the annotation group
            for k in self.anno_dict[anno]:
                segs.extend(self.get_segments(ch, k))
            # sort the segments
            segs.sort(key=lambda r: r[0])
            # check the sorting
            assert all(segs[i+1][0] >= segs[i][1] for i in range(len(segs)-1))
            # add lengths to total
            tot += sum(b-a for (a,b) in segs)
            # save to file
            f_save = self._savefile(anno, out_dir, ch)
            with open(f_save, 'w') as f:
                for (start, end) in segs:
                    line = '{} {} {}\n'.format(ch, start, end)
                    f.write(line)
        # print message for fraction of autosomes in anno
        print '{} {} {:.4f}'.format(self.cell_line, anno, tot/2.88e9)


def save_bed(savefile, chrom, segments):
    """save segments in a bed file format"""
    with open(savefile, 'w') as f:
        for (start, end) in segments:
            line = '{} {} {}\n'.format(chrom, start, end)
            f.write(line)

    return None


def merge_cell_lines(cell_lines, states):
    """merge states across 2 or more cell lines into composite annotation"""
    # load ChromHMMDatabase for each cell line once
    hmm_dict = {}
    for cl in cell_lines:
        hmm_dict[cl] = ChromHMMDatabase(cl)

    # for each state, merge data from each cell line
    for st in states:
        # create path to save merged files
        f_path = root_dir + '/data/bsanno/{}_merged'.format(st)
        if not os.path.isdir(f_path):
            os.mkdir(f_path)
        tot = 0
        # merge across cell lines one chromosome at a time
        for ch in human_autosomes:
            # create a master array for this chrom
            m_chrom = np.zeros(shape=chromosome_length(ch), dtype=bool)
            # get state data for each cell line
            for cl in cell_lines:
                # create empty chrom length array for annotations
                m = np.zeros(shape=chromosome_length(ch), dtype=bool)
                # get data for preloaded database
                hmm = hmm_dict[cl]
                # group data by annotation dictionary
                for k in hmm.anno_dict[st]:
                    # annotate each site in the empty chrom array
                    for (i,j) in hmm.get_segments(ch, k):
                        m[i:j] = 1
                # merge to the master array for the chrom
                m_chrom |= m

            # convert mask to segments
            segs = binary_mask_segments(m_chrom)

            # save segments
            f_save = f_path + '/{}.{}_merged.bed'.format(ch, st)
            save_bed(f_save, ch, segs)

            # keep a total of annotated sites across chroms
            tot += np.sum(m_chrom)

        # print the fraction of autosomal sites in state annotation
        print '{} merge: {}'.format(st, tot / 2.88e9)

    return None


def intersect_cell_lines(cell_lines, states):
    """merge states across 2 or more cell lines into composite annotation"""
    # load ChromHMMDatabase for each cell line once
    hmm_dict = {}
    for cl in cell_lines:
        hmm_dict[cl] = ChromHMMDatabase(cl)

    # for each state, intersect data from each cell line
    for st in states:
        # create path to save intersected files
        f_path = root_dir + '/data/bsanno/{}_intersected'.format(st)
        if not os.path.isdir(f_path):
            os.mkdir(f_path)
        tot = 0
        # intersect across cell lines one chromosome at a time
        for ch in human_autosomes:
            # create a master array for this chrom
            m_chrom = np.zeros(shape=chromosome_length(ch), dtype=bool)
            # get state data for each cell line
            for (i, cl) in enumerate(cell_lines):
                # create empty chrom length array for annotations
                m = np.zeros(shape=chromosome_length(ch), dtype=bool)
                # get data for preloaded database
                hmm = hmm_dict[cl]
                # group data by annotation dictionary
                for k in hmm.anno_dict[st]:
                    # annotate each site in the empty chrom array
                    for (i,j) in hmm.get_segments(ch, k):
                        m[i:j] = 1
                # use "OR" operator for first cell line
                if i == 0:
                    print 'i == 0'
                    m_chrom = m
                # use "AND" operator to intersect with the master array
                else:
                    m_chrom &= m

            # convert mask to segments
            segs = binary_mask_segments(m_chrom)

            # save segments
            f_save = f_path + '/{}.{}_intersected.bed'.format(ch, st)
            save_bed(f_save, ch, segs)

            # keep a total of annotated sites across chroms
            tot += np.sum(m_chrom)

        # print the fraction of autosomal sites in state annotation
        print '{} intersected: {}'.format(st, tot / 2.88e9)

    return None


def merge_peri_cds():
    """create composite of cds+periexonic segments"""
    p_open = root_dir + '/data/coords/nr'
    for ch in human_autosomes:
        f_p = p_open + '/peri/{}_knownGene_nonredundant_periSegments.bed'.format(ch)
        f_c = p_open + '/cds/{}_knownGene_nonredundant_cdsSegments.bed'.format(ch)
        segs = []
        for f_in in [f_p, f_c]:
            with open(f_in, 'r') as f:
                for line in f:
                    start, end = line.split()[1:3]
                    segs.append((int(start)-1, int(end)))
        segs.sort(key=lambda r: r[0])
        msegs = merge(segs)
        p_save = root_dir + '/data/bsanno/cds_peri_merge'
        if not os.path.isdir(p_save):
            os.mkdir(p_save)
        f_save = p_save + '/{}.cds_peri_merge.bed'.format(ch)
        save_bed(f_save, ch, msegs)

    return None


#%%

cell_lines = ['GM12878', 'H1hesc', 'Huvec', 'Hepg2', 'K562']
states = ['txn', 'pro', 'enh', 'ins', 'rep']
# merge_cell_lines(cell_lines, states)
intersect_cell_lines(cell_lines, states)

# merge_peri_cds()

# for cl in cell_lines:
#     hmm = ChromHMMDatabase(cl)
#     for st in states:
#         hmm.create_anno(st)

#%%