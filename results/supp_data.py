from classes.datastruct import DataStruct, np
from classes.mapstruct import MapStruct, _kwdpcalc, datetime
from classes.runstruct import root_dir
from functions import swap_root
from multiprocessing import Pool
from sys import argv, stderr, stdout
import os

__author__ = 'davidmurphy'


def main():
    # fdir = root_dir + '/result/final_files/sensitivity_tests'
    # f100 = fdir + '/YRI.pr95.conserved.BS1.6.CS0.0.170816142939.final.txt'
    # f1000 = fdir + '/YRI.pr95.conserved.1K.BS1.6.CS0.0.170816131908.final.txt'

    idir = root_dir + '/result/init_files'
    init = idir + '/YRI.pr95.cleanrun.BS1.6.CS0.0.NOT_STARTED.initial.txt'
    # init = idir + '/YRI.pr95.conserved.1K.BS1.6.CS0.0.170817043229.initial'
    # swap_root(init)
    d = DataStruct(init=init)
    for ch in d.chroms:
        d.chrom = ch
        s = d.substitution_counts(msk_file=d.neut_masks)
        fname = d.sub_count.replace('.npz', '.txt')
        header = 'upper_edge all_neutral neutral_substitutions'
        np.savetxt(fname=fname, X=s, fmt='%d %d %d', header=header)

    # m1 = MapStruct(init=init, pool=False)
    # m1.optimize()

    # p = Pool(processes=4)
    # pin = m1.pcalc_inputs(4)
    #
    # st = datetime.now()
    # message = 'starting {} -LLH calc... '.format(m1.process_type)
    # stderr.write(message)
    # stdout.flush()
    #
    # pout = p.map(_kwdpcalc, pin)
    # nlp = sum(pout) / m1.stat.pairs
    #
    # en = datetime.now()
    # message = '-LLH={}, time={}\n'.format(nlp, en - st)
    # stderr.write(message)
    # stdout.flush()
    #
    # # m1.log_likelihood(m1.params)
    #
    # m1.pool = None
    # m1.log_likelihood(m1.params)

    # m2 = MapStruct(init=f1000)
    # m2.vars.complete = False
    # fname = fdir + fstr
    # swap_root(fname)
    # m2 = MapStruct(init=fname)

    # m1.initial_stats()
    # m2.initial_stats()
    # m2.optimize()

    # params = np.array([-16.594535905414918631, -8.603346594444005291,
    #                    -8.628447304242518712, -8.301496479100499926,
    #                    -8.385001446528914215, -9.513574682513777248,
    #                    53.608687429280955428])
    # # params = m.params
    # m.vars.complete = False
    # res2 = m.log_likelihood(params=params)
    # print res2

    # m.initial_stats()
    # m.optimize()
    # m.initial_stats()
    # params = np.array([-16.594535905414918631, -8.603346594444005291,
    #                    -8.628447304242518712, -8.301496479100499926,
    #                    -8.385001446528914215, -9.513574682513777248,
    #                    -16.594535905414918631, -8.603346594444005291,
    #                    -8.628447304242518712, -8.301496479100499926,
    #                    -8.385001446528914215, -9.513574682513777248,
    #                    53.608687429280955428])
    # params = m.params
    # m.vars.complete = False
    # res1 = m.log_likelihood(params=params)
    # print res1

    # d = DataStruct(init=init_file, chrom='chr14')
    # d.compress_predictions()

    # cons = d.conservation_mask(wig_file=d.wigz_cons, percentile=0.80)
    # neut = d.conservation_mask(wig_file=d.wigz_ncons, percentile=0.0)
    # masked = neut & ~cons
    # print cons.sum(), neut.sum(), masked.sum()

    # mst = MapStruct(init=init_file)
    # mst.final_stats()
    # mst.save()
    # weighted tau match: 52.092399999999997817
    # unweighted tau match: 50.616810000000000969
    # wts = mst.nsites
    # for tau in np.arange(52.092, 52.095, 0.0001):
    #     mst.params[-1] = tau
    #     pred, pi0 = mst.calc()
    #     meanpred, meanpi0 = np.average(pred, weights=wts), np.average(pi0, weights=wts)
    #     print tau, mst.stat.meanpi-meanpred, 1-(mst.stat.meanpi/meanpi0)
    # mst.initial_stats()
    # dst = DataStruct(init=init_file, chrom='chr1')

    # for ch in dst.chroms[1:]:
    #     dst.chrom = ch
    #     con = np.loadtxt(dst.bs_targets.values()[0], usecols=(1, 2))
    #     seg = np.loadtxt(dst.genic_segments, usecols=(1, 2))
    #
    #     cons_mask = mask_segments(np.zeros(shape=dst.chlen, dtype=np.bool), con, flipoff=False)
    #     genic_mask = mask_segments(np.zeros(shape=dst.chlen, dtype=np.bool), seg, flipoff=False)
    #
    #     gcons = binary_mask_segments(mask=(cons_mask & genic_mask))
    #     with open(re.sub('Segments', 'genic', dst.bs_targets.values()[0]), 'w') as f:
    #         f.write('\n'.join('{}\t{}\t{}'.format(ch, i, j) for (i, j) in gcons))
    #
    #     ngcons = binary_mask_segments(mask=(cons_mask & ~genic_mask))
    #     with open(re.sub('Segments', 'nongenic', dst.bs_targets.values()[0]), 'w') as f:
    #         f.write('\n'.join('{}\t{}\t{}'.format(ch, i, j) for (i, j) in ngcons))


def cross_species_conserved(chrom, pct_cons, start=0):
    # initialize empty data struct
    d = DataStruct(make=True, chrom=chrom)

    # nodes are IDs for phylogenetic trees
    nodes = ['ape', 'primate', 'prosimian', 'euarchontoglires',
             'laurasiatheria', 'afrotheria', 'mammal', 'birds']
    scr = '/ifs/scratch/c2b2/gs_lab/dam2214/scratch'
    wigfile_template = '{}/wigz/{{s}}/{{s}}.{}.wig.gz'.format(scr, chrom)
    wigfiles = [wigfile_template.format(s='ape'), d.wigz_cons] + \
               [wigfile_template.format(s=s) for s in nodes[2:]]
    zipped = zip(nodes, wigfiles)[start:]

    # output template
    indent = max(map(len, nodes))
    stdout_template = '{{:<{m}}} {{:>10}} {{:>10}}\n'.format(m=indent)

    # generate a mask of non-conserved sites from 100-vertebrate alignment
    ncons = d.conservation_mask(fwig=d.wigz_ncons, pct=0)

    # print the header and the number of initial non-conserved sites
    header = '#spec_name tot_ncons spec_cons'.split()
    stderr.write(stdout_template.format(*header))
    stderr.write(stdout_template.format('init_ncons', ncons.sum(), 0))
    stdout.flush()

    for (n, f) in zipped:
        # reset cons ID corresponding to new file and generate top X% cons mask
        d.cons = n
        cons = d.conservation_mask(fwig=f, pct=pct_cons)
        # filter conserved sites from ncons
        ncons &= ~cons
        # print summaries for the current cons mask
        stderr.write(stdout_template.format(n, ncons.sum(), cons.sum()))
        stdout.flush()


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        # cross_species_conserved(chrom='chr22', pct_cons=0.80)
        # cons_ncons_filter(chrom='chr22', pct_cons='0.80')
        main()
    else:
        if len(argv) == 4:
            cross_species_conserved(argv[1], argv[2], int(argv[3]))
        else:
            print 'usage: cross_species_conserved <chrom> <pct_cons> <start>'
            exit()

