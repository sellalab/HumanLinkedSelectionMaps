__author__ = 'davidmurphy'


import os
import numpy as np
from itertools import izip
from sys import stderr, stdout, argv
from classes.knowngene import KnownGene
from classes.substitution import recursive_search
from classes.runstruct import ChromStruct, root_dir
from data_processing.data_tools import fasta_array, snpcount


def sample_substitutions(chrom, pop):
    """a function that probabilistically selects substitutions from SNP data"""
    # initialize a ChromStruct for file paths
    cst = ChromStruct(chrom=chrom, neut=pop)

    # load the HC ancestor for the chromosome - change '.' to 'N'
    anc = fasta_array(chrom, cst.ancs_files)
    anc[anc == '.'] = 'N'

    # load the reference for the chromosome
    ref = fasta_array(chrom, cst.refg_files)

    # print message for the number of ancestor only N sites
    anc_only_n = np.sum((anc == 'N') & (ref != 'N'))
    both_n = np.sum((anc == 'N') & (ref == 'N'))
    msg = 'anc_only_N={:.2e} both_N={:.2e}\n'.format(anc_only_n, both_n)
    stderr.write(msg)
    stdout.flush()

    # get initial nonmatching positions
    sub_pos_init = np.where(anc != ref)[0]

    # load the polymorphism data for the chromosome and population
    spos, rb, rc, ab, ac = snpcount(cst.snp_files, returnbases=True)
    # adjust the positions to be 0-based matching the arrays
    spos -= 1

    # randomly sample SNP alleles weighted by frequency (twice)
    n_alleles = rc+ac
    ref_prob = 1.0 * rc / n_alleles
    alt_prob = 1.0 * ac / n_alleles
    assert np.all(ref_prob+alt_prob == 1)
    pick_allele_1 = []
    pick_allele_2 = []
    for i in xrange(len(rb)):
        alleles = [rb[i], ab[i]]
        pick_1 = np.random.choice([0, 1], p=[ref_prob[i], alt_prob[i]])
        pick_2 = np.random.choice([0, 1], p=[ref_prob[i], alt_prob[i]])
        pick_allele_1.append(alleles[pick_1])
        pick_allele_2.append(alleles[pick_2])

    pick_allele_1 = np.array(pick_allele_1)
    pick_allele_2 = np.array(pick_allele_2)

    # identify subs by comparing modified ref to anc (first round)
    ref[spos] = pick_allele_1
    sub_pos_1 = np.where(anc != ref)[0]
    ref_1, anc_1 = ref[sub_pos_1], anc[sub_pos_1]
    out_1 = np.column_stack((sub_pos_1, ref_1, anc_1))

    # (second round)
    ref[spos] = pick_allele_2
    sub_pos_2 = np.where(anc != ref)[0]
    ref_2, anc_2 = ref[sub_pos_2], anc[sub_pos_2]
    out_2 = np.column_stack((sub_pos_2, ref_2, anc_2))

    # print some stats for the different samples of HC subs
    n = 1.0 * sub_pos_init.size
    p1 = np.sum(np.in1d(sub_pos_init, sub_pos_1)) / n
    p2 = np.sum(np.in1d(sub_pos_init, sub_pos_2)) / n
    msg = 'p1={:.4f} p2={:.4f}\n'.format(p1, p2)
    stderr.write(msg)
    stdout.flush()

    # save positions to compressed arrays
    f_token = '{}.{}.HC.subst.sample'.format(chrom, pop)
    f_out_1 = root_dir + '/data/div/hcder/{}.1.npz'.format(f_token)
    f_out_2 = root_dir + '/data/div/hcder/{}.2.npz'.format(f_token)
    np.savez_compressed(f_out_1, sub=out_1)
    np.savez_compressed(f_out_2, sub=out_2)

    return None


def get_annotations_from_sample(chrom, pop, sample, variant):
    """get the nonsynonymous substitutions from sampled substitutions"""
    # initialize a ChromStruct for file paths
    cst = ChromStruct(chrom=chrom, neut=pop)

    # load the list of substitutions
    f_token = '{}.{}.HC.subst.sample.{}'.format(chrom, pop, sample)
    f_in = root_dir + '/data/div/hcder/{}.npz'.format(f_token)
    sub = np.load(f_in)['sub']

    # load the reference for analyzing mRNA sequences
    href = fasta_array(chrom, cst.refg_files)
    # update the reference with the sampled alleles
    href[sub[:,0].astype(int)] = sub[:,1]
    # convert reference into a standard string for mRNA analysis
    href = ''.join(href)
    assert len(href) == cst.chlen

    # load the list of genes for the chromosome as KnownGene objects and cds
    search_list = []
    gene_list = {}
    with open(cst.files.fnrgene.format(ch=chrom), 'r') as f:
        f.next()
        for line in f:
            gene = KnownGene(line)
            entry = [(start, end, gene.name) for (start, end) in gene.cds]
            search_list.extend(entry)
            gene_list[gene.name] = gene

    # search for each substitution in the set of cds segments
    anno_list = []
    for (pos, ref, anc) in sub:
        # skip sites missing ancestor base
        if anc not in 'ACGT':
            continue
        pos = int(pos)
        search = recursive_search(pos, search_list, rtype='seg')
        # when a substitution is in a cds segment, classify the variant type
        if search is not None:
            gene = gene_list[search]
            variant_type = gene.variant(int(pos), anc, href)
            # if the variant matches variant parameter, add to list
            if variant_type == variant:
                anno_list.append(pos)
    anno_list = np.array(anno_list)

    if variant == 'NS':
        lbl = 'nonsyn'
    else:
        lbl = 'syn'

    # save file with all of the NS variants
    f_dir = root_dir + '/data/csanno/{}_{}_s{}'.format(pop, lbl, sample)
    if not os.path.isdir(f_dir):
        os.mkdir(f_dir)
    f_out = f_dir + '/{}.{}_{}_s{}.npz'.format(chrom, pop, lbl, sample)
    np.savez_compressed(f_out, pos=anno_list)

    return None


def main():
    # usage_1
    if len(argv) == 3:
        chrom, pop = argv[1:]
        sample_substitutions(chrom, pop)
    # usage_2
    elif len(argv) == 5:
        chrom, pop, sample, variant = argv[1:]
        get_annotations_from_sample(chrom, pop, sample, variant)
    else:
        print 'usage_1: sub_hc_substitutions <chrom> <pop>'
        print 'usage_2: sub_hc_substitutions <chrom> <pop> <sample> <variant>'
        exit(1)


if __name__ == '__main__':
    main()
