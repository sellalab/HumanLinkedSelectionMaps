__author__ = 'davidmurphy'

import os
import numpy as np
from sys import argv, stderr, stdout
from itertools import izip
from classes.runstruct import root_dir, ChromStruct
from likelihood.cllh_functions import predicted_pi
from precalc.lh_inputs import adjust_arrays


def n_alts(k, n):
    """unique het pairs with k alternates in sample size n"""
    return k * (n - k)


def n_pairs(n):
    """unique pairs for n alleles"""
    return 0.5 * n * (n-1)


def quad(a, b, c):
    """solutions to the quadratic formula"""
    s1 = (-b + np.sqrt((b**2 - 4*a*c))) / (2*a)
    s2 = (-b - np.sqrt((b**2 - 4*a*c))) / (2*a)
    return s1, s2


def c_term(npairs, pi):
    """c term for quadratic that solves alternative allele count from pi, n"""
    return pi * npairs


def simulated_data(pred, nsites, sample_size, random=False, verbose=True):
    """
    Return hom & het pair counts corresponding to predictions.

    When site counts per segment and sample size per site are fixed from the
    original data used to generate the predictions, this is essentially a
    deterministic simulation. The output for segment of the prediction map
    should be values of hom pairs and het pairs that yield the closest pi
    value to the prediction given the number of sites and samples per site.

    :param pred: predicted pi values in segments tiling the genome
    :param nsites: number of observed neutral sites in each segment
    :param sample_size: population sample size for neutral site data
    :param random: flag to use the more intensive random sampling routine
    :return: hom & het pair counts for each segment in pred
    :type pred: np.ndarray
    :type nsites: np.ndarray
    :type sample_size: int
    :type random: bool
    :rtype: np.ndarray
    """
    # set the number of possible pairs for each site based on sample size
    pairs = n_pairs(sample_size)

    # get the number of neutral segments from prediction map length
    nsegs = len(pred)

    # percentage completion printout to stderr
    pct = int(0.01 * nsegs)  # 1% of total grid points
    pct = max(pct, 1)
    if verbose:
        stderr.write('simulating neutral polymorphism data: ')
        stdout.flush()

    # simulate alleles in each neutral segment
    sim = np.zeros(shape=(nsegs, 2))
    for (idx, (cnt, prd)) in enumerate(izip(nsites.astype(int), pred)):
        # track percent progress
        if (idx % pct == 0) and verbose:
            msg = '.' if idx%(pct * 25) else '{}%'.format(idx/pct)
            stderr.write(msg)
            stdout.flush()

        # total allele pairs given neutral site counts
        tot_pairs = cnt * pairs

        # deterministic approximation: "solve" het, given  pi = het / tot_pairs
        if not random:
            # round to NEAREST int
            approx_het = int(prd * tot_pairs + 0.5)
            approx_hom = tot_pairs - approx_het
            sim[idx] = approx_hom, approx_het

        # simulate data with sampling noise
        else:
            # simple case for sample of 2 chromosomes
            if sample_size == 2:
                # mismatch pairs = sum of bernoulli random vars over sites
                approx_het = np.random.binomial(cnt, prd)
                approx_hom = tot_pairs - approx_het
                sim[idx] = approx_hom, approx_het

            # general case for n > 2
            else:
                # sample neutral freq spectrum for binomial SNP probability
                p_neut = np.random.beta(a=prd, b=prd, size=cnt)
                # p_neut = np.random.beta(a=prd, b=prd)

                # random sample of alt SNPS for each neutral site in block
                alts = np.random.binomial(sample_size, p_neut)
                # alts = np.random.binomial(sample_size, p_neut, cnt)

                # calc total hom/het pairs based on the number of alts
                s_hets = np.sum(n_alts(alts, sample_size))
                s_homs = tot_pairs - s_hets

                # record simulated hets/alts for the pred block
                sim[idx] = s_homs, s_hets

    # completed counter message and newline
    if verbose:
        stderr.write('\n')
        stdout.flush()

    return sim


def run_simulation(cst, sim_id, sample, params=None, random=True):
    """
    Use predictions to simulate poly data for set of params and save to file.

    :type cst: ChromStruct
    :param sim_id: suffix or name for output file
    :type sim_id: str
    :param params: parameters used to generate map of linked selection effects
    :type params: list
    """
    # load fixed arrays used to generate predicted map from params
    if cst.bnum:
        bs = np.load(cst.bs_files)['bs']
    else:
        bs = None
    if cst.cnum:
        cs = np.load(cst.cs_files)['cs']
    else:
        cs = None
    nt = np.load(cst.nt_files)['nt']
    dv = np.load(cst.dv_files)['dv']
    pl = np.load(cst.pl_files)['pl']
    nu = np.load(cst.nu_files)['nu']

    # # uniform mutation rate
    # nu = np.full(len(dv), cst.stat.meandiv)

    # load fixed arrays used to generate predicted map from params
    bs, cs, nu, _, dv, _, msk = adjust_arrays(cst, bs, cs, nu, nt, dv, pl)

    # # check that summaries and arrays make sense for the chromosomen
    # assert msk.size == cst.stat.totsegs[cst.cidx]
    # assert dv.size == cst.stat.msksegs[cst.cidx]

    # set params manually if new params provided otherwise use cst.params
    if params is not None:
        cst.params = params

    # get a map of predictions for the given params and precalc maps
    pred_map = predicted_pi(cst.params, cst, nu, bs, cs)

    # # simulate polymorphism data for each segment with neutral sites
    # poly_sim = simulated_data(pred_map, dv, cst.stat.indv)
    poly_sim = simulated_data(pred_map, dv, sample, random)

    # create an empty array for all segments tiling the chrom
    sim = np.zeros(shape=(msk.size, 2))

    # use mask to sort simulated polymorphic data into correct segments
    sim[msk] = poly_sim

    # get simulated data file path and save sim data sorted into chroms
    fsim = cst.nt_simfiles(sim_id).format(ch=cst.chrom)
    np.savez_compressed(fsim, nt=sim)

    return None


if __name__ == '__main__':
    # local
    if root_dir.startswith('/Users/davidmurphy'):
        fdir = '{}/result/final_files/sims'.format(root_dir)
        fname = '00-YRI.pr95.clustinf.initial.BS1.6.CS0.0.180408214759.final.txt'
        finit = '{}/{}'.format(fdir, fname)
        chst = ChromStruct('chr21',  init=finit)

        # from classes.reformat_files import reformat_fields
        # reformat_fields(finit)
        # chst = ChromStruct('chr21', init=finit)
        # chst.files.__init__(chst.dict)
        chst.stat.indv = 216
        sid = 'infprm.dtrm'
        prm = chst.params

        run_simulation(chst, sid, prm)

    else:
        if len(argv) != 5:
            print 'usage: neutpoly_sim <ch> <init> <sim_id> <params>'
            exit(1)

        ch, init, sid = argv[1:4]
        prm = eval(argv[4])
        chst = ChromStruct(chrom=ch, init=init)
        # patch for working with older files since in the new cst
        # chst.files.__init__(chst.dict)

        run_simulation(chst, sid, prm)
