__author__ = "davidmurphy"

# TODO: finish member functions for this class

chr_length = dict(chr1=249250621, chr2=243199373, chr3=198022430,
                  chr4=191154276, chr5=180915260, chr6=171115067,
                  chr7=159138663, chr8=146364022, chr9=141213431,
                  chr10=135534747, chr11=135006516, chr12=133851895,
                  chr13=115169878, chr14=107349540, chr15=102531392,
                  chr16=90354753, chr17=81195210, chr18=78077248,
                  chr19=59128983, chr20=63025520, chr21=48129895,
                  chr22=51304566, chrY=59373566, chrX=155270560)

genetic_code = dict(TAG='-', TAA='-', TGA='-',
                    GCA='A', GCC='A', GCG='A', GCT='A',
                    TGT='C', TGC='C',
                    GAC='D', GAT='D',
                    GAG='E', GAA='E',
                    TTT='F', TTC='F',
                    GGT='G', GGG='G', GGA='G', GGC='G',
                    CAT='H', CAC='H',
                    ATC='I', ATA='I', ATT='I',
                    AAG='K', AAA='K',
                    CTT='L', CTG='L', CTA='L', CTC='L', TTA='L', TTG='L',
                    ATG='M',
                    AAC='N', AAT='N',
                    CCT='P', CCG='P', CCA='P', CCC='P',
                    CAA='Q', CAG='Q',
                    AGG='R', AGA='R', CGA='R', CGC='R', CGG='R', CGT='R',
                    AGC='S', AGT='S', TCT='S', TCG='S', TCC='S', TCA='S',
                    ACC='T', ACA='T', ACG='T', ACT='T',
                    GTA='V', GTC='V', GTG='V', GTT='V',
                    TGG='W',
                    TAT='Y', TAC='Y')


def complement_base(b):
    # mini-function to complement a base
    complement_code = dict(A='T', C='G', G='C', T='A')
    return complement_code[b]


def complement_strand(s):
    # mini-function to complement a string of bases
    return ''.join(complement_base(b) for b in s)


def recursive_search(pos, segments):
    """
    Mini-function to find the a segment containing a site
    :param pos: query site
    :param segments: a set of segments
    :return segment[m]: the segment index containing the site (if any)
    """
    n = len(segments)
    m = n / 2
    if n == 0:
        return None
    if segments[m][0] <= pos <= segments[m][1]:
        return segments[m]
    elif pos < segments[m][0]:
        return recursive_search(pos, segments[:m])
    elif pos > segments[m][1]:
        return recursive_search(pos, segments[m + 1:])


class KnownGene(object):
    """
    An object that represents a UCSC known gene.
    ============================================

    The UCSC Genes track is a set of gene predictions based on programs from RefSeq, GenBank, CCDS, Rfam, and the tRNA 
    Genes track. The track includes both protein-coding genes and non-coding RNA genes. Both types of genes can produce
    non-coding transcripts, but non-coding RNA genes do not produce protein-coding transcripts. This is a moderately
    conservative set of predictions. Transcripts of protein-coding genes require the support of one RefSeq RNA, or
    one GenBank RNA sequence plus at least one additional line of evidence. Transcripts of non-coding RNA genes
    require the support of one Rfam or tRNA prediction. Compared to RefSeq, this gene set has generally about 10%
    more protein-coding genes, approximately four times as many putative non-coding genes, and about twice as many
    splice variants.


    Data fields:
    ============

    field           example             description
    -----           -------             -----------
    name            uc001aaa.3	        Name of gene
    chrom           chr1	 	        Reference sequence chromosome or scaffold
    strand 	        +	                + or - for strand
    txStart 	    11873	 	        Transcription start position
    txEnd 	        14409	            Transcription end position
    cdsStart 	    11873               Coding region start
    cdsEnd 	        11873               Coding region end
    exonCount 	    3                   Number of exons
    exonStarts 	    11873,12612,13220,  Exon start positions
    exonEnds 	    12227,12721,14409,	longblob 	  	Exon end positions
    proteinID       BOB11               UniProt display ID, UniProt accession, or RefSeq protein ID
    alignID 	    uc001aaa.3          Unique identifier (GENCODE transcript ID for GENCODE Basic

    (source: http://genome.ucsc.edu/cgi-bin/hgTables)

    """

    def __init__(self, datastring):
        """
        The class is initialized from a line in the knownGene.txt file
        :param datastring: '\t' separated values for the knownGene
        """

        # keep the programs string for writing to file later on
        self.datastring = datastring

        # strip "\n", split on "\t", read in each field as a string
        (name,
         chrom,
         strand,
         txStart,
         txEnd,
         cdsStart,
         cdsEnd,
         exonCount,
         exonStarts,
         exonEnds,
         proteinID,
         alignID) = datastring.strip().split("\t")

        # some can just remain strings
        self.name, self.chrom, self.strand = name, chrom, strand
        self.proteinID, self.alignID = proteinID, alignID

        # check that strand and chrom make sense (for human chromosomes only)
        assert self.strand in ["+", "-"]
        # assert self.chrom in ["chr{}".format(c) for c in range(1, 23) + ["X", "Y", "M"]] --> not true, some odd ones

        # other variables are converted to int, etc.
        (self.txStart,
         self.txEnd,
         self.cdsStart,
         self.cdsEnd,
         self.exonCount) = map(int, [txStart, txEnd, cdsStart, cdsEnd, exonCount])

        # starts and ends are split by ","
        exonStarts = map(int, exonStarts.split(",")[:-1])
        exonEnds = map(int, exonEnds.split(",")[:-1])

        # some constant rules that must be true
        assert all(exonStarts[i] < exonEnds[i] for i in xrange(self.exonCount))
        assert (self.cdsStart >= exonStarts[0] and self.cdsEnd <= exonEnds[-1])

        # turn the starts and ends into a list of segment (start, stop) tuples
        self.exonSegments = [(exonStarts[i], exonEnds[i]) for i in xrange(self.exonCount)]
        self.exonLength = sum(s[1] - s[0] for s in self.exonSegments)

        # get the introns from the inner segments between exons
        self.intronSegments = [(exonEnds[i], exonStarts[i + 1]) for i in xrange(self.exonCount - 1)]
        self.intronLength = sum(s[1] - s[0] for s in self.intronSegments)

        assert self.intronLength + self.exonLength == self.txEnd - self.txStart

        # some boolean tests of the gene:
        self.coding = cdsEnd != cdsStart

        # get the coding segments if they exist as well as any UTRs
        if self.coding:

            self.utrSegments = []
            self.utr5 = []
            self.utr3 = []

            # find the exon segments where cds starts and ends
            startExon = recursive_search(self.cdsStart, self.exonSegments)
            endExon = recursive_search(self.cdsEnd, self.exonSegments)

            # if they are in one exon, go no further
            if startExon == endExon:
                self.codingSegments = [(self.cdsStart, self.cdsEnd)]

            else:
                # use this info to create 5' and 3' cds ends
                cds5prime = (self.cdsStart, startExon[1])
                cds3prime = (endExon[0], self.cdsEnd)

                # indices of middle 5' and 3' exons
                m5pr = self.exonSegments.index(startExon) + 1
                m3pr = self.exonSegments.index(endExon)

                # string together the coding segments
                self.codingSegments = [cds5prime] + self.exonSegments[m5pr:m3pr] + [cds3prime]

            # some other info about the gene
            self.codingCount = len(self.codingSegments)
            self.codingLength = sum(s[1] - s[0] for s in self.codingSegments)
            self.translated = bool(proteinID)
            self.truncated = bool(self.codingLength % 3)

            # get the 5' UTR for this gene
            if self.txStart < self.cdsStart:
                # find the index of the last segment before the cds start
                utr5end = self.exonSegments.index(startExon)
                utr5 = self.exonSegments[0:utr5end] + [(startExon[0], self.cdsStart)]

                # only keep segments which have a non-zero length
                self.utr5 = [(x[0], x[1]) for x in utr5 if x[1] > x[0]]

            # get the 3' UTR for this gene
            if self.txEnd > self.cdsEnd:
                # find index of the first segment beyond the cds end
                utr3start = self.exonSegments.index(endExon) + 1
                utr3 = [(self.cdsEnd, endExon[1])] + self.exonSegments[utr3start:]

                # only keep segments which have a non-zero length
                self.utr3 = [(x[0], x[1]) for x in utr3 if x[1] > x[0]]

            # add the segments to the list of UTRs in increasing order, depending on strand
            if self.strand == "+":
                self.utrSegments = self.utr5 + self.utr3
            else:
                # if "-" strand, swap everything
                tmp3 = self.utr5
                self.utr5 = self.utr3
                self.utr3 = tmp3
                self.utrSegments = self.utr3 + self.utr5

            self.utrLength = sum(s[1] - s[0] for s in self.utrSegments)
            self.utrCount = len(self.utrSegments)

            # utrs don't overlap, all segments sum to gene length, at most utr + cds has 2 more segments than exons
            assert all(self.utrSegments[i][1] <= self.utrSegments[i + 1][0] for i in xrange(self.utrCount - 1))
            assert self.exonLength == self.codingLength + self.utrLength
            assert self.exonCount <= self.codingCount + self.utrCount <= self.exonCount + 2

            # get 500bp of TSS
            if self.strand == '+':
                self.tss = [(self.txStart - 500, self.txStart)]
                # based on the strand, tss will be either the lowest or greatest position in the gene
                self.geneStart = self.tss[0][0]
                self.geneEnd = self.txEnd
            else:
                self.tss = [(self.txEnd, self.txEnd + 500)]
                self.geneStart = self.txStart
                self.geneEnd = self.tss[0][1]

            # get the peri-exonic regions (100bp up and downstream of exons)
            self.splicing = []
            if self.exonCount > 1:
                for i in xrange(self.exonCount - 1):
                    # if the spacing is > 200 create two separate segments
                    if self.intronSegments[i][1] - self.intronSegments[i][0] > 200:
                        lower = self.intronSegments[i][0], self.intronSegments[i][0] + 100
                        upper = self.intronSegments[i][1] - 100, self.intronSegments[i][1]
                        self.splicing.append(lower)
                        self.splicing.append(upper)
                    # if the spacing is less than or equals 200 just take the entire intron
                    else:
                        self.splicing.append(self.intronSegments[i])
                # assert non-overlapping and sorted
                assert all(self.splicing[i][1] <= self.splicing[i + 1][0] for i in xrange(len(self.splicing) - 1))

    @property
    def exons(self):
        """return a list of exon segments from the current gene"""
        return self.exonSegments

    @property
    def introns(self):
        """return a list of intron segments from the current gene"""
        return self.intronSegments

    @property
    def cds(self):
        """return a list of CDS (coding) segments from the current gene"""
        return self.codingSegments

    # @property
    # def utr5(self):
    #     """return the UTR5 segment from the current gene"""
    #     return self.utr5
    #
    # @property
    # def utr3(self):
    #     """return the UTR3 segment from the current gene"""
    #     return self.utr3

    # def __str__(self):
    #     """String representation of the knownGene"""
    #     out = "{:<10}\t{:<10}\n{:<10}\t{:<10}\n{:<10}\t{:<10}"
    #     return out.format("Chromosome", self.chrom, "Coding", str(self.coding), "Gene", self.proteinID)

    def __len__(self, part="mRNA"):
        """
        Get the length of the transcript or mRNA
        :param part: specify whether to return length of the mRNA or transcript
        """
        if part == "mRNA":
            return self.codingLength
        elif part == "transcript":
            return self.exonLength

    def mrna(self, ref):
        """
        Get the mRNA sequence for the gene
        :param ref: reference genome as a string
        :return mrna: the sequences of bases (with U = T) from the ref genome
        """
        if self.truncated:
            print "Warning: gene {} has a truncated sequence".format(self.name)

        mrna = ""
        for (start, stop) in self.codingSegments:
            mrna += ref[start:stop]

        return mrna

    def transcript(self, ref):
        """
        Get the transcript sequence for the gene
        :param ref: reference genome as a string
        :return transcript: the sequences of bases (with U = T) from the ref genome
        """
        transcript = ""
        for (start, stop) in self.exonSegments:
            transcript += ref[start:stop]

        return transcript

    def protein(self, ref):
        """
        Translate the mrna sequence to a protein sequence
        :param ref: reference genome as a string
        :return protein: single letter amino acid coded protein
        """
        protein = ""
        mrna = self.mrna(ref).upper()  # up-case for safety
        if self.strand == "+":
            for i in xrange(0, self.codingLength, 3):
                codon = mrna[i:i + 3]
                protein += genetic_code[codon]

        elif self.strand == "-":
            # reverse and complement the mrna for '-' strand
            mrna = "".join(map(complement_base, mrna[::-1]))
            for i in xrange(0, self.codingLength, 3):
                codon = mrna[i:i + 3]
                protein += genetic_code[codon]

        return protein

    def variant(self, pos, alt, ref, aa=False):
        """
        Determine the effects of a SNP on the protein.
        NOTE: subtract 1 if input is from a 1-based coordinate system
        :param pos: position where the SNP occurred --> NOTE: convert between 0 vs. 1 based coords
        :param alt: the new base at pos to compare with reference mrna
        :param ref: the reference genome for the chrom where the gene is found
        :param aa: optional flag to return the translation effect (i.e. Q>V, etc.)
        :return effect: the effect of the SNP --> S, NS, STOP
        :return translation: if aa=True, translate the change caused by the SNP
        """

        # get the mrna for the gene
        seq = self.mrna(ref)

        # find the segment and segment index where the SNP occurs
        seq_segment = recursive_search(pos, self.codingSegments)
        segment_idx = self.codingSegments.index(seq_segment)

        # calculate the position in the mrna sequence where the SNP occurs
        if segment_idx == 0:
            # if it is in the first segment just subtract lower range from pos
            segment_pos = pos - self.codingSegments[0][0]
        else:
            # if it is in a segment other than 0, add the number of bases per segment before with the number of bases
            # into the segment where pos is found
            segment_pos = sum([self.codingSegments[i][1] - self.codingSegments[i][0] for i in xrange(segment_idx)])
            segment_pos += pos - self.codingSegments[segment_idx][0]

        # insert the alternate base to the new sequence
        new_seq = seq[:segment_pos] + alt + seq[segment_pos + 1:]

        # sequences must not change length
        assert len(seq) == len(new_seq)

        # check the new codon
        codon_start = (segment_pos / 3) * 3
        old_codon = seq[codon_start:codon_start + 3].upper()  # up-case for safety
        new_codon = new_seq[codon_start:codon_start + 3].upper()

        # if len(old_codon) != 3:
        #     print self.datastring
        #     raise KeyError(pos, alt)

        if self.strand == "+":
            old_aa = genetic_code[old_codon]
            new_aa = genetic_code[new_codon]
        else:
            # reverse and complement for '-' strand
            old_codon = "".join(map(complement_base, old_codon[::-1]))
            new_codon = "".join(map(complement_base, new_codon[::-1]))
            old_aa = genetic_code[old_codon]
            new_aa = genetic_code[new_codon]

        if new_aa == "-" and new_aa != old_aa:
            return "STOP"
        elif new_aa != old_aa:
            return "NS"
        else:
            return "S"
