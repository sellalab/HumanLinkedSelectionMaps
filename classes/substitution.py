import numpy as np


__author__ = "davidmurphy"


def recursive_search(pos, segments, rtype=None):
    """
    Mini-function to search for a position in a set of coding segments
    :param pos: query site
    :param segments: a set of segments
    :param rtype: control whether the segment where the position was found or simply 'True' if within some segment in
    the search set
    :return segment[m]: the segment index containing the site (if any)
    """
    # TODO: can we truncate the bottom of the list as we go to make it shorter?

    n = len(segments)
    m = n / 2
    if n == 0:
        return None
    if segments[m][0] <= pos < segments[m][1]:  # <-- subs must be LESS than the upper range of the segment
        if rtype == 'pos':
            return True  # for looking in arbitrary segments
        elif rtype == 'seg':
            return segments[m][2]  # for gene annotation (the 3rd position in the segment is some hash or index)
        else:
            return segments[m]
    elif pos < segments[m][0]:
        return recursive_search(pos, segments[:m], rtype=rtype)
    elif pos > segments[m][1]:
        return recursive_search(pos, segments[m + 1:], rtype=rtype)


class Substitution(object):
    """
    This class represents a single SNP and its genomic annotation
    """

    def __init__(self, datastring):

        # save the programs string to write out later
        self.datastring = datastring.strip()

        # split up the programs to init member variables
        (self.chrom,
         pos,
         self.derived,
         self.ancestor) = datastring.strip().split("\t")

        # initialize member variable "pos" as an int
        self.pos = int(pos)

    def __str__(self):
        return self.datastring

    def annotate(self, tagged_segments, genedict, ref):
        """
        This function uses a reference genome of KnownGene objects stored in a dictionary to determine if the
        substitution occurs in a gene. If it occurs in a gene, it is analyzed to determine whether it is S, NS or STOP
        :param tagged_segments: a dict with chromosomes as keys and lists of all KnownGene codingSegments with a
        gene.name tag for a given chromosome
        """
        # the tag is the gene.name that is returned if the substitution is in a gene
        tag = recursive_search(self.pos, tagged_segments)

        if tag:
            # if a tag is returned, get the gene out of the genedict
            gene = genedict[tag]
            self.classification = gene.variant(self.pos, self.ancestor, ref)
        else:
            # otherwise the substitution is not in any genes so it is labelled noncoding
            self.classification = 'noncoding'


def inorout(positions, segments):
    """
    Determine which positions fall within a set of segments
    :param positions: Genomic 1-based positions
    :param segments: 1-based segments
    :return:
    """

    insegments = []
    for p in positions:
        if recursive_search(pos=p, segments=segments):
            insegments.append(p)

    return insegments


# debug
# from gzip import open as zopen
# c = 22
# ddir = '/Users/davidmurphy/GoogleDrive/linked_selection/data'
#
# dsubs_file = '{}/subs/derived_substitutions/chr{}_derived_substitutions_hg19.coords.gz'.format(ddir, c)
# with zopen(dsubs_file, 'r') as f:
#     dsubs = [int(line.split()[1]) for line in f if not line.startswith('#')]
#
# cons_file = '{}/pyLS/bsanno/primate_cons95_Segments/chr{}.primate_cons95_Segments.bed'.format(ddir, c)
# segs = np.loadtxt(cons_file, usecols=(1, 2))
#
# inp1 = inorout(positions=dsubs, segments=segs)
# print 'chr{} complete.'.format(c)

# outf = '{}/pyLS/subs/chr{}_primate_cons95_substitutions'.format(data, c)
# pos = np.loadtxt('{}/subs/annovar/chr{}.subst.variant_function'.format(data, c), usecols=(3, ), dtype=int)
# segs1 = np.loadtxt('{}/coords/nr/cons/primate/segs/chr{}_primate_cons95_Segments.bed'.format(data, c),
#                    usecols=(1, 2), dtype=int)
# segs2 = np.loadtxt('/Volumes/Macintosh HD/Users/davidmurphy/GoogleDrive/linked_selection/data/coords/nr/cds'
#                    '/chr{}_knownGene_nonredundant_cdsSegments.bed'.format(c), usecols=(1, 2),
#                    dtype=int)

# GETTING ALL THE SUBSTITUTIONS IN PRIMATE CONSERVED

# subs = np.array(inp1)
# np.save(outf, subs)

# inp2 = inorout(positions=pos, segments=segs2)
# print 'chr{}: cons={}, exons={}'.format(c, len(inp2), '1')

