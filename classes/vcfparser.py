from gzip import open as gzopen

__author__ = "davidmurphy"


class VCFReader(object):
    """
    A class for processing VCF records using by parsing data with helpful subclasses
    """
    def __init__(self, vcf_file):

        # attributes initialized null or empty
        self._line = ''
        self._fileheader = []

        # create an open file from the vcf file, init the closed flag to False
        self._vcf = gzopen(vcf_file, 'r')
        self._closed = False

        # read the full header from the file
        while not self.line.startswith('#CHROM'):
            self._getline()
            self._fileheader.append(self.line)

    def __iter__(self, external_ids=None):
        """yield Record objects for data lines until the file is closed"""
        while not self.closed:
            yield self.getrecord()

    def _getline(self):
        """update the current line in the vcf file and then return it"""
        try:
            self._line = self._vcf.next()
        except StopIteration:
            self.close()  # when the last line is reached, close the file
            raise StopIteration

    def getrecord(self):
        """return Record object initialized with the current line"""
        self._getline()
        return Record(self.line)

    def close(self):
        self._vcf.close()
        self._closed = True

    @property
    def ids(self):
        return self.header[9:]

    @property
    def line(self):
        return self._line

    @property
    def closed(self):
        return self._closed

    @property
    def header(self):
        """just the HEADER line of the VCF file (i.e. the line that begins with #CHROM) split into columns"""
        return self._fileheader[-1].split('\t')

    @property
    def metainfo(self):
        """all of the double commented meta-information lines of the VCF file"""
        return self._fileheader[:-1]

    @property
    def fileheader(self):
        """fileheader lines as '\n'-separated string"""
        return ''.join(self._fileheader)


class Record(object):
    """
    A light-weight mapping of VCF data elements to class properties. 
    """
    def __init__(self, line):
        """split line of VCF data into elements"""
        self.data = line[:-1].split('\t')  # take up to index len(data)-1 to drop '\n' char
        assert len(self.data) >= 9

    def __str__(self):
        return '\t'.join(self.data)

    @property
    def sample(self):
        return 2 * (len(self.data) - 9)

    @property
    def chrom(self):
        return self.data[0]

    @property
    def pos(self):
        return int(self.data[1])

    @property
    def vid(self):
        return self.data[2]

    @property
    def ref(self):
        return self.data[3]

    @property
    def alt(self):
        return self.data[4]

    @property
    def qual(self):
        return self.data[5]

    @property
    def filter(self):
        return self.data[6]

    @property
    def info(self):
        """return a dict of ID=VALUE pairs for all info IDs in the current record"""
        return dict(tuple(x.split('=')) if '=' in x else (x, True) for x in self.data[7].split(';'))

    @property
    def vformat(self):
        """return a list of format IDs for the current record"""
        return self.data[8].split(':')

    @property
    def idata(self):
        """return list of dict(format_id=value) for each individual in the current record"""
        if len(self.data) > 9:
            return [dict((k, v) for (k, v) in zip(self.vformat, indv.split(':'))) for indv in self.data[9:]]

    @property
    def allele_counts(self):
        """return a the reference and alternate counts for the record, if applicable"""
        alts = sum(int(d['GT'][0]) + int(d['GT'][-1]) for d in self.idata)  # sum the 1s across genotypes
        refs = self.sample - alts  # total alleles - alts = refs
        return refs, alts
