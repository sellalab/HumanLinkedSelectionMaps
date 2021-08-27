from knowngene import KnownGene, chr_length
from collections import defaultdict
from itertools import ifilter
from nrgenes import nrgenes
from merge import merge
import re


__autor__ = "davidmurphy"


# directories
coords_dir = "/Users/davidmurphy/GoogleDrive/linked_selection/data/coords/"


def mask(segment):
    """
    Mini-function to use with builtin filter function, use to remove pathological length 1 segments
    :param segment: a START, END tuple representing a segment
    :return boolean: evaluates True for segments greater than length 1, else false
    """
    return segment[1] - segment[0] > 1


def gene_data(ucsc_file):
    """
    This function writes 5 files from programs in the UCSC knownGene annotations based on criteria in function nrgenes.
    See Hernandez et al. 2011 supplement for details on these criteria. This implementation allows exons to overlap
    as long as coding sequences do not overlap.

    Files written by this function:
    -------------------------------

    (1) a set of genes in ucsc knownGene format for non-overlapping genes
    (2) a set of positions for the midpoint of each codingSegment in the set for a given chromosome
    (3) a set of exonSegment coordinates in bed format
    (4) a set of codingSegment coordinates in bed format
    (5) a set of utrSegment coordinates in bed format
    (6) a set of intronSegment coordinates in bed format
    (7) a log file with stats about the run
    (8) a master log for stats of the whole set of autosomes

    :param ucsc_file: the path for the UCSC knownGene file to be used
    """

    genes = defaultdict(list)
    # TODO: this is just temporarily changed to get X programs on its own
    # autosomes = ["chr{}".format(c) for c in xrange(1, 23)]
    chromosomes = ['chrX']

    # some global stats for the run
    auto_bases        = 0
    auto_genes        = 0
    auto_exons        = 0
    auto_introns      = 0
    auto_cds          = 0
    auto_utrs         = 0
    auto_exon_bases   = 0
    auto_coding_bases = 0
    auto_utr_bases    = 0
    auto_intron_bases = 0

    with open(ucsc_file, "r") as f:
        f.next()  # skip header
        for line in f:
            g = KnownGene(line)
            # keep all genes that are coding and non-truncated
            if g.coding and not g.truncated:
                # store each set of genes for a given chromosome as a list in dict
                genes[g.chrom] += [g]

    # process each gene list
    for ch in chromosomes:

        # take each list of genes and process the list with gene_set to get a non-redundant list
        nr_genes  = nrgenes(genes[ch])

        midpoints = []  # get all midpoints for all codingSegments
        exons     = []  # get all exonSegments
        cds       = []  # get all codingSegments
        utrs      = []  # get all UTR 5' & 3'
        introns   = []  # get all introns

        for gene in nr_genes:

            # calculate the midpoint
            pts = [s[0] + (s[1] - s[0]) / 2 for s in gene.codingSegments]
            midpoints += pts

            # gather intron, exon, UTR and coding segments, mask out pathological length 1 segments
            exons   += filter(mask, gene.exonSegments)
            cds     += filter(mask, gene.codingSegments)
            utrs    += filter(mask, gene.utrSegments)
            introns += filter(mask, gene.intronSegments)

        # NOTE: Where (start + 1 == end), the segment is simply dropped from the stats and the output bed file
        exonic_bases = sum(s[1] - s[0] for s in exons)
        coding_bases = sum(s[1] - s[0] for s in cds)
        utr_bases    = sum(s[1] - s[0] for s in utrs)
        intron_bases = sum(s[1] - s[0] for s in introns)

        # fraction to write out:
        exonic_base_fraction  = "{:.3f}".format(1.0 * exonic_bases / chr_length[ch])
        coding_bases_fraction = "{:.3f}".format(1.0 * coding_bases / chr_length[ch])
        utr_bases_fraction    = "{:.3f}".format(1.0 * utr_bases / chr_length[ch])
        intron_bases_fraction = "{:.3f}".format(1.0 * intron_bases / chr_length[ch])

        # check that ordering is correct
        assert all(exons[i][1] <= exons[i + 1][0] for i in xrange(len(exons) - 1))
        assert all(cds[i][1] <= cds[i + 1][0] for i in xrange(len(cds) - 1))
        assert all(utrs[i][1] <= utrs[i + 1][0] for i in xrange(len(utrs) - 1))
        assert all(introns[i][1] <= introns[i + 1][0] for i in xrange(len(introns) - 1))

        # increment global counters
        auto_bases        += chr_length[ch]
        auto_genes        += len(nr_genes)
        auto_exons        += len(exons)
        auto_introns      += len(introns)
        auto_cds          += len(cds)
        auto_utrs         += len(utrs)
        auto_exon_bases   += exonic_bases
        auto_coding_bases += coding_bases
        auto_utr_bases    += utr_bases
        auto_intron_bases += intron_bases

        # output (1)
        f1 = "{path}nr/gene/{chrom}_knownGene_coding_nonredundant_genes.txt".format(path=coords_dir, chrom=ch)
        with open(f1, "w") as o1:
            o1.write("# {} UCSC knownGene coding nonredundant genes\n".format(ch))
            o1.write("".join(g.datastring for g in nr_genes))  # note: datastring contains "\n" already

        # output (2)
        f2 = "{path}nr/pts/{chrom}_knownGene_nonredundant_coding_midpoints.txt".format(path=coords_dir, chrom=ch)
        with open(f2, "w") as o2:
            o2.write("# {} UCSC coding segment midpoints from nonredundant coding gene set\n".format(ch))
            o2.write("\n".join(str(pt) for pt in midpoints))

        # ** NOTE **
        # ----------
        # to conform with "BED" style coordinates, the 0 - based segment (start, end) is converted to (start + 1, end)
        # in the following output files. Where (start + 1 == end), the segment is simply dropped.

        # output (3)
        f3 = "{path}nr/exon/{chrom}_knownGene_nonredundant_exonSegments.bed".format(path=coords_dir, chrom=ch)
        with open(f3, "w") as o3:
            o3.write("# {} UCSC exon segments from nonredundant coding gene set\n".format(ch))
            o3.write("\n".join("{}\t{}\t{}".format(ch, start + 1, stop) for (start, stop) in exons))

        # output (4)
        f4 = "{path}nr/cds/{chrom}_knownGene_nonredundant_codingSegments.bed".format(path=coords_dir, chrom=ch)
        with open(f4, "w") as o4:
            o4.write("# {} UCSC coding segments from nonredundant coding gene set\n".format(ch))
            o4.write("\n".join("{}\t{}\t{}".format(ch, start + 1, stop) for (start, stop) in cds))

        # output (5)
        f5 = "{path}nr/utr/{chrom}_knownGene_nonredundant_utrSegments.bed".format(path=coords_dir, chrom=ch)
        with open(f5, "w") as o5:
            o5.write("# {} UCSC UTR 5' and 3' segments from nonredundant coding gene set\n".format(ch))
            o5.write("\n".join("{}\t{}\t{}".format(ch, start + 1, stop) for (start, stop) in utrs))

        # output (6)
        f6 = "{path}nr/intron/{chrom}_knownGene_nonredundant_intronSegments.bed".format(path=coords_dir, chrom=ch)
        with open(f6, "w") as o6:
            o6.write("# {} UCSC intron segments from nonredundant coding gene set\n".format(ch))
            o6.write("\n".join("{}\t{}\t{}".format(ch, start + 1, stop) for (start, stop) in introns))

        # output (7)
        f7 = "{path}nr/stats/{chrom}_knownGene_nonredundant_stats.txt".format(path=coords_dir, chrom=ch)
        with open(f7, "w") as o7:
            out1 = "{:>18} {:>11,}\n"
            out2 = "{:>18} {:>11}\n"
            header    = "Statistics for {} UCSC knownGene coding nonredundant genes:".format(ch)
            o7.write("{}\n{}\n".format(header, "-" * len(header)))
            o7.write(out1.format("Genes:", len(nr_genes)))
            o7.write(out1.format("Exons:", len(exons)))
            o7.write(out1.format("Introns:", len(introns)))
            o7.write(out1.format("Coding segments:", len(cds)))
            o7.write(out1.format("UTR segments:", len(utrs)))
            o7.write(out1.format("Autosomal bases:", chr_length[ch]))
            o7.write(out1.format("Exonic bases:", exonic_bases))
            o7.write(out1.format("Intronic bases:", intron_bases))
            o7.write(out1.format("Coding bases:", coding_bases))
            o7.write(out1.format("UTR bases:", utr_bases))
            o7.write(out2.format("Fraction exonic:", exonic_base_fraction))
            o7.write(out2.format("Fraction intronic:", intron_bases_fraction))
            o7.write(out2.format("Fraction coding:", coding_bases_fraction))
            o7.write(out2.format("Fraction UTR:", utr_bases_fraction))

    # # output (8)
    # f8 = "{path}stats/knownGene_nonredundant_autosome_stats.txt".format(path=coords_dir)
    # with open(f8, "w") as o8:
    #     out1 = "{:>18} {:>13,}\n"
    #     out2 = "{:>18} {:>13}\n"
    #     header    = "Statistics for UCSC knownGene coding nonredundant genes on autosomes:"
    #     o8.write("{}\n{}\n".format(header, "-" * len(header)))
    #     o8.write(out1.format("Genes:", auto_genes))
    #     o8.write(out1.format("Exons:", auto_exons))
    #     o8.write(out1.format("Introns:", auto_introns))
    #     o8.write(out1.format("Coding segments:", auto_cds))
    #     o8.write(out1.format("UTRs:", auto_utrs))
    #     o8.write(out1.format("Autosomal bases:", auto_bases))
    #     o8.write(out1.format("Exonic bases:", auto_exon_bases))
    #     o8.write(out1.format("Intronic bases:", auto_intron_bases))
    #     o8.write(out1.format("Coding bases:", auto_coding_bases))
    #     o8.write(out1.format("UTR bases:", auto_utr_bases))
    #     o8.write(out2.format("Fraction exonic:", "{:.3f}".format(1.0 * auto_exon_bases / auto_bases)))
    #     o8.write(out2.format("Fraction intronic:", "{:.3f}".format(1.0 * auto_intron_bases / auto_bases)))
    #     o8.write(out2.format("Fraction coding:", "{:.3f}".format(1.0 * auto_coding_bases / auto_bases)))
    #     o8.write(out2.format("Fraction UTR:", "{:.3f}".format(1.0 * auto_utr_bases / auto_bases)))

    return None


def merged_exons(ucsc_file, truncated=False, coding=False):
    """
    This function returns a merged set of non-overlapping exon ranges for each chromosome in the ucsc_file,
    with possible filters applied AND the set of genes used to make them. This contrasts with the above function that
    only keeps a SINGLE set of genes that are coding, non-truncated and non-overlapping (without any merges),
    with annotated start and stop codons, no internal stops, where overlaps trigger a choice of the longest coding
    sequence in the overlapping genes.
    :param ucsc_file: knownGene.txt file (required due to class usage)
    :param truncated: boolean flag for whether or not to filter out truncated genes
    :param coding: boolean flag for whether or not to filter out non-coding annotated genes.
    :return: (1) a set of ranges written to bed file, some of which have been modified due to merging. because of the merge,
    these are not "real" exons but rather exonic regions (2) the genes used to construct the merged exons
    """

    # make a list of exonic region for each autosome stored in a dict
    genes = defaultdict(list)
    autosomes = ["chr{}".format(c) for c in xrange(1, 23)]

    with open(ucsc_file, "r") as f:
        f.next()  # skip header
        for line in f:
            g = KnownGene(line)

            # all filters
            if coding is True and truncated is True and g.coding is True and g.truncated is False:
                genes[g.chrom].append(g)

            # coding filter only
            elif coding is True and truncated is False and g.coding is True:
                genes[g.chrom].append(g)

            # no filters
            elif coding is False and truncated is False:
                genes[g.chrom].append(g)

            # error: cannot be truncated if not coding
            elif coding is False and truncated is True:
                print "error: gene cannot be truncated if it is not coding"
                exit(1)

    # now go through each autosome and merge the set of segments for each gene on the autosome
    for ch in autosomes:

        # list to store all exon segments
        exons = []

        # iterate through all genes
        for gene in genes[ch]:

            # collect all exon segments into one list
            exons += gene.exonSegments

        # sort all of the segments in order from start position
        exons.sort(key=lambda s: s[0])

        # merged returns a generator, filter it with ifilter to remove length 1 segments
        merge_exons = ifilter(mask, merge(exons))

        # check sum of segment spans
        # print "{}: {}".format(ch, sum(b - a for a, b in merge_exons))

        # save all the genes that were used to construct the merged exon segments
        f1 = "{path}merged/gene/{chrom}_knownGene_coding_merged_genes.txt".format(path=coords_dir, chrom=ch)
        with open(f1, "w") as o1:
            o1.write("# {} UCSC knownGene coding merged genes\n".format(ch))
            o1.write("".join(g.datastring for g in genes[ch]))

        f2 = "{path}merged/exon/{chrom}_knownGene_coding_merged_exonSegments.bed".format(path=coords_dir, chrom=ch)
        with open(f2, "w") as o2:
            o2.write("# {} UCSC merged exon segments\n".format(ch))
            o2.write("# Filters: coding-only = {}, non-truncated = {}\n".format(coding, truncated))
            # NOTE: Where (start + 1 == end), the segment is simply dropped.
            o2.write("\n".join("{}\t{}\t{}".format(ch, a + 1, b) for (a, b) in merge_exons))

    return None


def near_exons(chromosome, geneset):
    """
    Return segments 500bp upstream of the TSS and 100bp 5' & 3' of exon start & end positions (based on
    recommendations of Jonathan Pritchard)
    :param chromosome: human chromosome (use to create filename)
    :param geneset: flag for the dataset to use (currently nr or merged)
    :return: write the segments out to a bed file
    """

    assert geneset in ["nr", "merged"]
    ch = "chr{}".format(chromosome)

    # final list to write to file
    near_exon_segments = []

    # find the ucsc file for that chromosome (nonredundant gene set at the moment)
    ucsc_file = "{path}{gset}/gene/{ch}_knownGene_coding_nonredundant_genes.txt"
    ucsc_file = ucsc_file.format(path=coords_dir, ch=ch,gset=geneset)

    # convert genes into KnownGene objects as programs list
    with open(ucsc_file, "r") as f:
        _genes = (KnownGene(x) for x in f if not x.startswith("#"))

        # take the first gene to check back if the 500bp 5' of the next gene overlaps
        gene1 = next(_genes)

        # get the 500bp 5' of first TSS
        near_exon = [(max(1, gene1.exonSegments[0][0] - 500), gene1.exonSegments[0][0])]

        # get the peri-exonic regions for this gene
        for (start, end) in gene1.intronSegments:
                # take up to 100bp downstream of the segment or as much as possible until next exon
                near_exon.append((start, min(start + 100, end)))
                # take up to 100bp upstream of the segment or as much as possible between exons
                near_exon.append((max(end - 100, start), end))
        assert len(near_exon) <= 2 * gene1.exonCount - 1

        # save the tail position in order to check that the 500bp 5' of TSS does not overlap last gene
        tail = gene1.exonSegments[-1][1]
        near_exon_segments.extend(near_exon)

        # for each gene, get intronic segments up to 100bp upstream of exon starts and 100 bp downstream of exon ends
        for gene in _genes:
            near_exon = [(max(tail, gene.exonSegments[0][0] - 500), gene.exonSegments[0][0])]
            for (start, end) in gene.intronSegments:
                near_exon.append((start, min(start + 100, end)))
                near_exon.append((max(end - 100, start), end))

            assert len(near_exon) <= 2 * gene.exonCount - 1
            tail = gene.exonSegments[-1][1]
            near_exon_segments.extend(near_exon)

    # quick fix for overlaps due to the -500bp from TSS step -- sort, remove len == 1 segments, then merge
    # TODO: maybe fix this or reconsider how it should be done at some point?
    near_exon_segments.sort(key=lambda a: a[0])
    near_exon_segments = filter(mask, near_exon_segments)
    new_segments = merge(near_exon_segments)

    # write to a new bed file
    f1 = "{path}{gset}/near/{ch}_knownGene_nonredundant_nearExonSegments.bed"
    f1 = f1.format(path=coords_dir, ch=ch, gset=geneset)
    with open(f1, "w") as o1:
        o1.write("# {} UCSC merged gene set near-exon segments\n".format(ch))
        o1.write("# 500bp upstream of TSS, 100bp intronic bases 5' and 3' of each exon segment\n")
        # NOTE: Where (start + 1 == end), the segment is simply dropped.
        o1.write("\n".join("{}\t{}\t{}".format(ch, a + 1, b) for (a, b) in new_segments))

    return None


def genic_segments(ucsc_file, outdir):
    """
    Pull down tss, utr5, cds, splice, utr3 segments from coding genes and merge together for a single set of "genic"
    segments
    :param ucsc_file: knownGene file
    :param outdir: path to send new files
    """
    segments = defaultdict(list)
    with open(ucsc_file, "r") as f:
        _ = f.readline()  # skip header
        for line in f:
            g = KnownGene(line)
            # keep coding autosomes only, index lists of segments by chrom
            if g.coding and not g.truncated and g.chrom == 'chrX':
                # the components of the "genic" annotation (exon = utr5/3 + cds)
                segments[g.chrom].extend(g.tss + g.exonSegments + g.splicing)

    # write programs to file one chrom at a time
    for ch in segments:
        # turn segment lists into merged-segment iterators and save to new file
        merged = merge(segments=sorted(segments[ch], key=lambda a: a[0]))
        newfile = '{d}/{c}.genicMerge.bed'.format(d=outdir, c=ch)
        with open(newfile, 'w') as o:
            o.write('# Genic segments: tss + exons + splicing, merged and sorted -- chr, start, end\n')
            o.write('\n'.join('{}\t{}\t{}'.format(ch, start + 1, end) for start, end in ifilter(mask, merged)))

    return None


def merge_conserved(file_names):
    """
    A program that is used to merge together two or more BED format set of annotations (e.g., making one "gcons" from
    McVicker's gcons_ex and gcons_nex). The new annotations are returned after merging
    :param file_names: a list of BED files to be merged -- these must have the same chromosome indicated in the name
    :type file_names: list
    :return merged: a generator that yields sorted segments from the merge result
    :rtype merged: generator
    """
    # combine all of the segments from all of the files in one list
    segments = []
    for afile in file_names:
        f = open(afile, 'r')
        for line in f:
            start, end = map(int, line.split()[1:])
            segments.append((start - 1, end))  # subtract 1 to set back to 0-based coordinates
        f.close()
    # sort segments by segment-start position
    segments.sort(key=lambda t: t[0])
    # return a merged iterator
    return merge(segments=segments)


def main():
    ucsc = '/Users/davidmurphy/GoogleDrive/linked_selection/data/ch_features/knownGene.txt'
    gdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/coords/genicMerge'
    # gcons = '/Users/davidmurphy/GoogleDrive/linked_selection/data/coords/gcons'
    #
    # # merging McVicker conserved sites
    # for ch in xrange(2, 23):
    #     input_files = [gcons + s for s in '/ex/chr{}_gcons_ex.bed'.format(ch), '/nex/chr{}_gcons_nex.bed'.format(ch)]
    #     output_file = gcons + '/segs/chr{}_gcons_segments.bed'.format(ch)
    #     joined = merge_conserved(file_names=input_files)
    #     with open(output_file, 'w') as f:
    #         f.write('\n'.join('chr{}\t{}\t{}'.format(ch, start + 1, end) for start, end in joined) + '\n')
    #     print '{:5>} complete'.format('chr' + str(ch))

    genic_segments(ucsc_file=ucsc, outdir=gdir)
    # merged_exons(ucsc, coding=True)
    # for c in xrange(1, 2):
    #     near_exons(chromosome=c, geneset="nr")

    # run:
    # from datetime import datetime as dt
    # st = dt.now()
    # merged_exons(ucsc, coding=True)
    # gene_data(ucsc)
    # print "time = {}".format(dt.now() - st)

    # from datetime import datetime as dt
    # st = dt.now()
    # for-loop1 time=0:00:01.276823
    # for-loop2 time=0:00:01.189824
    # list-comp time=0:00:01.332121

if __name__ == '__main__':
    main()
