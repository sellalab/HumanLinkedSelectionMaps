from sys import stdin
from subprocess import call
from gzip import open as gzopen
from classes.dictlike import DictLike

__author__ = 'davidmurphy'


class AlignmentSequence:
    """
    A class representing "s" lines (i.e. sequence lines) in the MAF file format
    """
    def __init__(self, line):
        # save a copy of the data string
        self._str = line

        # split on single whitespace
        data = line.split()
        assert len(data) == 7

        # initialize members
        self.src = data[1]
        self.start = data[2]
        self.size = data[3]
        self.strand = data[4]
        self.srcSize = data[5]
        self.text = data[6]

        # get just the species tag from src
        # self.species = data[1].split('.')[0]

    def __str__(self):
        return self._str

    def __len__(self):
        return len(self.text)

    def setstring(self):
        """update the internal data_string with the current attribute values"""
        self._str = 's '+' '.join([self.src, self.start, self.size,
                                   self.strand, self.srcSize, self.text])+'\n'

    @property
    def species(self):
        return self.src.split('.')[0]

    @property
    def chrom(self):
        return self.src.split('.')[1]


class Parser:
    """
    A class to manage maf lines, file names and species sets, keep track of the N*L bases in the data and start new
    files and processes when N*L crosses some threshold
    """
    def __init__(self, species_string):
        """
        Initialize a new writer with the set of species to keep from MAF sequence lines 
        :param species_string: a string of ',' separated species
        """
        # split individual species into a set
        self._species = set(species_string.split(','))

        # current header and sequences
        self._current_block = []

    def __len__(self):
        """return current number of sequences in block (don't count header)"""
        return len(self._current_block) - 1

    def get_block(self, maf_lines, as_string=True):
        """
        A generator that yields processed blocks of MAF data according to
        preset selectors. The generator processes MAF lines from an open
        file or stdin and yields processed blocks when block-separator empty
        '\n' lines are read
        :param maf_lines: iterable source of MAF format lines
        :param as_string: flag to return string block vs. list of block data
        """
        for line in maf_lines:
            # header line: start new block beginning with the header
            if line.startswith('a'):
                self._current_block = [line]

            # sequence line: store AlignmentSequence if ID in self.species
            elif line.startswith('s'):
                species_id = line.split()[1].split('.')[0]
                if species_id in self._species:
                    new_sequence = AlignmentSequence(line)
                    self._current_block.append(new_sequence)

            # block-separator: yield current block
            elif line.startswith('\n'):
                # yield the current block if >= nonref sequence
                if len(self) > 1:
                    # default output: string
                    if as_string:
                        yield self.block_string
                    else:
                        yield self.block_data

            # check that other lines are in 'ieq' categories like they should be
            else:
                assert line[0] in '#ieq'

    def _maskhg19(self):
        """replace reference seq with following seq. keep reference coords"""
        if len(self._current_block) > 2:
            self._current_block[0].text = self._current_block[1].text
            self._current_block[0].size = self._current_block[1].size
            self._current_block[0].setstring()
            self._current_block.remove(self._current_block[1])
        else:
            self._current_block = []

    @property
    def block_data(self):
        """return header and AlignmentSequence for each seq in the block"""
        return self._current_block

    @property
    def block_string(self):
        """return block as string: header and seqs separated by newlines"""
        return ''.join(str(b) for b in self._current_block)

    @property
    def block_header(self):
        """header line associated with MAF block (first position of list)"""
        return self._current_block[0]


class Splitter:
    """
    An class for splitting large MAF files into smaller subunits and removing unused lines. The split_file function
    is an iterator which writes files up to the specified sequence legnth and then yields the name of the just
    completed MAF file.
    """
    def __init__(self, chrom, taxon, species, out_dir, char_limit=None,
                 maf_file=None):
        """
        Load input and ouput paths
        :param chrom: chromosome for the current data
        :param taxon: taxon name given to splitter output
        :param species: comma-separated string of species IDs used by Parser to select MAF lines
        :param out_dir: destination for MAF files 
        :param char_limit: file size limit in characters (i.e. bytes). Default is None (i.e., parse the whole file)
        :param maf_file: path to MAF file - if no maf_file is given, Splitter reads from stdin
        """
        self.chrom = chrom
        self.taxon = taxon
        self.parser = Parser(species_string=species)
        self.out_dir = out_dir
        self.char_limit = char_limit
        self.maf_file = maf_file

        # counters for chars and files
        self.char_count = 0
        self.file_num = 1

        # current open write file
        self.current_file = open(self.current_filename, 'w')
        # MAF data source: if no file is given, read from stdin
        self.maf_lines = gzopen(maf_file, 'r') if self.maf_file else stdin

    def split_file(self):
        """
        A generator that splits a full chromosome MAF file into sub-files. Each time a file is completed the 
        generator yields the name of the completed file.
        """
        # process lines into blocks with Parser until EOF triggers StopIteration
        while self.maf_lines:
            try:
                # rest counters and open new file at the top of the loop AFTER
                # the most recent yield
                if self._stop:
                    self._yield(new_file=True)
                # try to get next block from Parser and write to current file
                block_string = self.parser.get_block(self.maf_lines).next()
                self.current_file.write(block_string)
                # update char count for the current file
                self.char_count += len(block_string)
                # if char count crosses limit, yield current file name start new file
                if self._stop:
                    yield self.current_filename

            except StopIteration:
                self._yield(new_file=False)
                yield self.current_filename

    def _yield(self, new_file):
        """Close current file. If new_file=True: update counters and open a new file; new_file=False: close all"""
        self.current_file.close()
        if new_file:
            self.char_count = 0
            self.file_num += 1
            self.current_file = open(self.current_filename, 'w')
        else:
            self.maf_lines.close()
            self.maf_lines = None

    @property
    def _stop(self):
        """return True if the Splitter should stop and create a new file; always False if no char_limit"""
        if self.char_limit is not None:
            return self.char_count > self.char_limit
        else:
            return False

    @property
    def current_filename(self):
        """Return write file name reflecting the current file number"""
        if self.char_limit is not None:
            return '{}/{}.{}.n{}.maf'.format(self.out_dir, self.chrom,
                                             self.taxon, self.file_num)
        else:
            return '{}/{}.{}.maf'.format(self.out_dir, self.chrom, self.taxon)


class PhastCaller(DictLike):
    """
    A class which splits large MAF files or MAF stdin streams into small files for processing in phastCons
    """
    def __init__(self, chrom, model, char_limit, info_string, maf_file=None, **kwargs):
        super(PhastCaller, self).__init__()

        # phastcons root and essential file paths
        # self.root = '/Users/davidmurphy/GoogleDrive/linked_selection/myprograms/cluster_code/phast/'
        root = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch/phast/'
        self.logs = root + 'logs'
        self.mafs = root + 'mafs'
        self.wigs = root + 'wigs'
        self.jobs = root + 'jobs'

        # phastCons exe path
        self.prog = '/ifs/data/c2b2/gs_lab/dam2214/phastcons/phast_13/bin/phastCons'

        # required command line inputs
        self.chrom = chrom
        self.fnum = 1
        self.model = model
        self.taxon, species = info_string.split(':')

        # splitter will yield new files to send to phastCons
        self.splitter = Splitter(self.chrom, self.taxon, species, char_limit, self.mafs, maf_file)
        self.maf = self.splitter.current_filename
        self.log = '{logs}/{chrom}.{taxon}.n{{fnum}}.log'.format(**self.dict)
        self.job = '{jobs}/{chrom}.{taxon}.n{{fnum}}.sh'.format(**self.dict)
        self.wig = '{wigs}/{chrom}.{taxon}.n{{fnum}}.wig'.format(**self.dict)

        # default tuning params
        self.rho = 0.3
        self.ecov = 0.3
        self.elen = 45

    def send_jobs(self, call_limit=None):
        """
        Iterate through 'call_limit' MAF files with the splitter and analyze with phastcons
        :param call_limit: number of calls to make to phastCons
        """
        for f in self.splitter.split_file():
            # if using a call limit, break once call limit is reached
            if (call_limit is not None) and (self.fnum > call_limit):
                break
            self.maf = f
            self._jobfile()
            # print self.qsub_cmmd
            call(self.qsub_cmmd, shell=True)
            self.fnum += 1

    def _jobfile(self):
        """write a temp shell file to run via qsub"""
        job = self.job.format(fnum=self.fnum)
        with open(job, 'w') as f:
            f.write('#!/bin/sh\n' + self.phast_cmmd + self.cleanup_cmmd)

    @property
    def phast_cmmd(self):
        """return phastCons command line arguments"""
        temp = '{prog} -R {rho} -C {ecov} -E {elen} -N {chrom} -i MAF {maf} {model} > {wig}\n'.format(**self.dict)
        return temp.format(fnum=self.fnum)

    @property
    def cleanup_cmmd(self):
        """returns rm command to delete MAF file after it has been used by phastCons"""
        return 'rm {}\n'.format(self.maf)

    @property
    def qsub_cmmd(self):
        """return the qsub shell command"""
        temp = 'qsub -l mem=1G,time=:5: -cwd -j y -o {log} {job}'.format(**self.dict)
        return temp.format(fnum=self.fnum)
