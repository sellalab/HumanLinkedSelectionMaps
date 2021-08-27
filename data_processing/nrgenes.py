import gzip
from collections import defaultdict
from classes.knowngene import KnownGene
from classes.runstruct import ChromStruct

__author__ = "davidmurphy"


def gene_dict():
    """return dict of lists of knownGene genes hashed by chrom"""

    # get ucsc file from chromstuct and process each line with KnownGene class
    genes = defaultdict(list)
    cst = ChromStruct(chrom='chr1')
    ucsc_file = cst.files.fgene
    chroms = list(cst.chroms) + ['chrX']
    with open(ucsc_file, "r") as f:
        f.next()  # skip header
        for line in f:
            g = KnownGene(line)
            # keep all genes that are coding, non-truncated and map to X or Aut
            if g.coding and (not g.truncated) and (g.chrom in chroms):
                # build lists of genes for each of the chroms
                genes[g.chrom].append(g)

    return genes


def nrgenes(genelist):
    """
    A function to create a list of non-redundant genes. When multiple genes
    overlap, the largest of the overlapping genes is kept. Truncated genes
    (i.e. codingLength % 3 != 0) are removed. Genes with internal stop codons
    or missing start or stop are also removed.
    See Hernandez et al., 2011 supplemental info: "Genic and CNC annotations"
    :param genelist: a list of KnownGene objects
    :return nr: the non-redundant list of genes
    """

    # check that all chroms agree
    ch = genelist[0].chrom
    assert all([g.chrom == ch for g in genelist])

    # get genome file from chromstruct
    cst = ChromStruct(chrom=ch)

    # truncated and noncoding should not even be put into the list
    assert not any([g.truncated for g in genelist])
    assert all([g.coding for g in genelist])

    # load the reference genome
    genome_file = cst.refg_files
    with gzip.open(genome_file, "r") as f:
        f.next()  # skip header
        sequence = f.read().replace("\n", "")

    # remove genes that do not translate correctly
    good_genes = []
    for gene in genelist:

        # translate to protein
        prot = gene.protein(sequence)
        # internal stop
        if "-" in prot[:-1]:
            continue
        # missing stop
        elif prot[-1] != "-":
            continue
        # missing start
        elif prot[0] != "M":
            continue
        else:
            good_genes += [gene]

    # ** NOTE **
    # change this (01/14/16) so that sorting is done with respect to exonStart
    # sort the list by the first position of the first codingSegment
    # good_genes.sort(key=lambda g: g.codingSegments[0][0])
    good_genes.sort(key=lambda g: g.exonSegments[0][0])

    # put the first gene in to initialize the list
    nr = [good_genes[0]]
    for gene in good_genes[1:]:

        # the current last gene in the list
        current = nr[-1]

        # all segments must be in order (overlap is possible for fully contiguous exons in the cds so <= is used)
        assert all([gene.codingSegments[i][1] <= gene.codingSegments[i + 1][0] for i in xrange(gene.codingCount - 1)])

        # if the gene overlaps with the last gene, replace shorter gene with longer gene
        # ** NOTE **
        # change this (01/14/16) so that overlap is called for EXON overlap, not CDS overlap
        # change this (05/29/16) so that TSS may not overlap previous gene
        # if gene.codingSegments[0][0] < current.codingSegments[-1][1]:
        if gene.geneStart <= current.geneEnd:
            if gene.codingLength > current.codingLength:
                nr[-1] = gene

        # if it does not overlap, append it to the list
        else:
            nr += [gene]

    # ** NOTE **
    # change this (01/14/16) so that overlap check is with respect to exonStart
    # assert all([nr[i].codingSegments[-1][1] < nr[i+1].codingSegments[0][0] for i in xrange(len(nr)-1)])
    assert all([nr[i].exonSegments[-1][1] <= nr[i + 1].exonSegments[0][0] for i in xrange(len(nr) - 1)])

    return nr


def nr_data(anno):
    """
    Get a dictionary of ordered non-redundant exons for each chrom.
    :param anno: annotation to get info from (e.g., exon)
    exon start
    :return nrex: a list of non-redundant exon coordinates
    """

    pass




def main():
    cst = ChromStruct(chrom='chr1')
    gdict = gene_dict()
    # save nr exons to files
    for ch in gdict.keys():
        # reset chromosome in chromstruct for correct file names
        cst.chrom = ch
        # get the non-redundant gene list
        nr = nrgenes(gdict[ch])
        # write non-redundant exons
        with open(cst.nr_exons, 'w') as f:
            # write the header
            hdr = '# {} UCSC exon segments from nonredundant coding gene set'
            f.write(hdr.format(ch))
            # for each gene, write each line of exon coordinates in BED format
            for g in nr:
                fmt = '{}\t{}\t{}'
                # set generator of exon segment strings for the gene
                ex_str = (fmt.format(ch, i+1, j) for (i, j) in g.exonSegments)
                # write exon segment strings separated by newlines
                f.write('\n' + '\n'.join(s for s in ex_str))


if __name__ == '__main__':
    main()