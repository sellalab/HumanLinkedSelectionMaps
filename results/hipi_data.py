from sys import argv

import matplotlib.pyplot as plt

from classes.hipidata import MapStruct, HiPiData, np, os
from functions import swap_root

# from merge import merge

__author__ = 'davidmurphy'


def main_remote():
    init1 = argv[1]
    init2 = argv[2]
    snps = []

    for init_file in [init1, init2]:
        swap_root(init_file, other_strings=[('margs_', 'op'), ('self.', 'initializer.')])
        hdd = HiPiData(init=init_file)
        snps.append(hdd.snps)

    matching = np.in1d(snps[0], snps[1])
    fraction_match = np.sum(matching, dtype='f8') / (0.5 * sum(map(len, snps)))
    print 'fraction of matching SNPs = {}'.format(fraction_match)


def main_local():
    # files
    data_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/data'
    res_dir = data_dir + '/pyLS/result/'
    # init = res_dir + 'final_files/filter.only.keep.BGC.BS1CS0.170406065527.final.txt'
    # init = res_dir + 'final_files/prim-95-cons-div.BS1CS0.170406062432.final.txt'
    gaps = dict(chr10=[39254935, 42254935],
                chr11=[51644205, 54644205],
                chr12=[34856694, 37856694],
                chr13=[16000000, 19000000],
                chr14=[16000000, 19000000],
                chr15=[17000000, 20000000],
                chr16=[35335801, 38335801],
                chr17=[22263006, 25263006],
                chr18=[15460898, 18460898],
                chr19=[24681782, 27681782],
                chr1=[121535434, 124535434],
                chr2=[92326171, 95326171],
                chr3=[90504854, 93504854],
                chr4=[49660117, 52660117],
                chr5=[46405641, 49405641],
                chr6=[58830166, 61830166],
                chr7=[58054331, 61054331],
                chr8=[43838887, 46838887],
                chr9=[47367679, 50367679],
                chr20=[26369569, 29369569],
                chr21=[11288129, 14288129],
                chr22=[13000000, 16000000])
    # for k in gaps:
    #     gaps[k].sort(key=lambda t: t[0])
    #     gaps[k] = list(merge(gaps[k]))
    #     lst = []
    #     for i in xrange(len(gaps[k]) - 1):
    #         lst.extend([gaps[k][i][1], gaps[k][i + 1][0]])
    #     gaps[k] = lst

    tags = [' non-BGC', '']
    alphas = [1.0, 0.6]
    inits = [res_dir + 'final_files/filter.all.BGC-2.BS1CS0.170404005648.final.txt',
             res_dir + 'final_files/filter.only.keep.BGC.BS1CS0.170406065527.final.txt']
    i = 0
    plt.figure(figsize=(13, 7))
    for init in inits[1:]:
        hdd = HiPiData(init=init, seg=True)
        # snps = []
        # for init_file in [init1, init2]:
        #     hdd = HiPiData(init=init_file)
        #     snps.append(hdd.snps)
        # m1 = m2 = 0
        # n = 0
        # n1 = n2 = 0
        # for (ch1, ch2) in zip(snps[0], snps[1]):
        #     matching1 = np.in1d(ch1[:, 0], ch2[:, 0])
        #     matching2 = np.in1d(ch2[:, 0], ch1[:, 0])
        #     # print len(ch1), len(ch2)
        #     n1 += len(ch1)
        #     n2 += len(ch2)
        #     m1 += np.sum(matching1)
        #     m2 += np.sum(matching2)
        #     n += len(ch2)
        # fraction_match1 = 1.0 * m1 / n1
        # fraction_match2 = 1.0 * m2 / n2
        # print n1, n2
        # print 'match={}, tot={}, fraction={}'.format(m1, n1, fraction_match1)
        # print 'match={}, tot={}, fraction={}'.format(m2, n2, fraction_match2)
        # exit(1)

        # # distribution of high diversity across chroms

        # all sites
        pos = hdd.pos.astype('f8')
        # rng = hdd.ranges.astype('f8')
        sites = hdd.nsites
        ch = hdd.chrom_tag

        # normalize relative positions on chroms
        for c in xrange(1, 23):
            chrom = 'chr{}'.format(c)
            cidx = np.where(ch == c)[0]
            cpos = pos[cidx]
            cmax, cmin = cpos.max(), cpos.min()
            pos[cidx] = 100.0 * (cpos - cmin) / (cmax - cmin)
            icen = gaps[chrom][0]
            jcen = gaps[chrom][1]
            ipos = np.where(cpos < icen)[0]
            jpos = np.where(cpos > jcen)[0]
            cpos[ipos] -= cmin
            cpos[ipos] /= (icen - cmin)
            cpos[jpos] -= cmax
            cpos[jpos] /= (jcen - cmax)
            # pos[cidx] = (100.0 * cpos)

            # print '{}\t0\t{}'.format(chrom, rng[cidx].min() + 5e6)
            # print '{}\t{}\t{}'.format(chrom, rng[cidx].max() - 5e6, chromosome_length(chrom))
            # p = pos[cidx]  # / (float(chromosome_length(chrom)) * 0.01)
            # s = sites[cidx]
            # high pi sites
            # hdd_cidx = np.where(ch[hdd.hidx] == c)[0]
            # hdd_pos = pos[hdd.hidx][hdd_cidx]  # / (float(chromosome_length(chrom)) * 0.01)
            # hdd_sites = sites[hdd.hidx][hdd_cidx]

        hdd_pos = pos[hdd.hidx]
        hdd_sites = sites[hdd.hidx]

        plt.hist(pos, bins=500, weights=sites, label='all data' + tags[i], alpha=1, normed=1)
        plt.hist(hdd_pos, bins=500, weights=hdd_sites, label='high diversity data' + tags[i], alpha=0.9,
                 normed=1)
        plt.xlabel('relative chromosomal position', fontsize=16, labelpad=10)
        plt.xticks(fontsize=16)
        plt.ylabel('site fraction', fontsize=16, labelpad=10)
        plt.yticks(fontsize=16)
        i += 1

    plt.legend()
    plt.show()

    # plt.figure(figsize=(13, 7))
    # est = MutationRatesEstimator(make=True)
    # err = est.estimator_relative_error(newmin=1e3)
    # plt.hist(err, bins=100, range=(0, 1), alpha=0.75, label='new {}kb'.format(int(est.fixed.muest_window / 1e3)))
    # print np.mean(err), np.median(err)
    #
    # est.err_rate = 'tmp'
    # est.sub_count = [cnt.format(ch) for ch in est.chroms]
    # err = est.estimator_relative_error(1e3)
    # plt.hist(err, bins=100, range=(0, 1), alpha=0.75, label='old {}kb'.format(int(est.fixed.muest_window / 1e3)))
    # print np.mean(err), np.median(err)
    #
    # plt.legend()
    # plt.show()
    # exit(1)

    # est = MutationRatesEstimator(make=True)
    # plt.figure(figsize=(13, 7))
    # bases = np.concatenate([est.substitution_counts(ch) for ch in est.chroms])[:, 1]
    # bases = bases[bases > 0]
    # plt.hist(bases, bins=1000, label='new', normed=1, alpha=0.75)
    # print np.mean(bases), np.median(bases)
    # bases = np.concatenate([np.load(cnt.format(ch)) for ch in est.chroms])[:, 1]
    # bases = bases[bases > 0]
    # plt.hist(bases, bins=1000, label='old', normed=1, alpha=0.75)
    # print np.mean(bases), np.median(bases)
    # plt.legend()
    # plt.show()

    # # analyze some mapping quality features
    # labdict = dict(ps='called bases', phcnt='multispecies alignment', acgt='ACGT', het='het pairs')
    # # features = []
    # features = ['het']
    # for feat in features:
    #     label = '% {} in segment'.format(labdict[feat])
    #     hddhist, gnhist = hdd.compare_feature(feature=feat, remove_self=1, random_sample=0)
    #     distribution_plot(hddhist, hdd.nsites, gnhist, hdd.nsites, label)

    # # analyze repeat composition for major repeat types
    # repeat_types = 'Low_complexity Simple_repeat DNA LTR LINE SINE'.split()
    # repeat_types = ['SINE']
    # repeat_types = ['filtered']
    # for repeat in repeat_types:
    #     label = '% {} in segment'.format(repeat.lower().replace('_', ' '))
    #     hddhist, gnhist = hdd.compare_repeat(mst=mst, repeat=repeat, random_sample=0)
    #     distribution_plot(hddhist, hdd.nsites, gnhist, mst.nsites, label)

    return None


# noinspection PyTypeChecker
def distribution_plot(hddhist, hddsites, gnhist, gnsites, label):
    """
    A plot of the ratio of % high diversity / genome neutral sites in bins reflecting some segment feature
    :param hddhist: histogram counts for high diversity data
    :param hddsites: vector of neutral sites per high diversity segment
    :param gnhist: histogram counts for genome data
    :param gnsites: vector of neutral sites per genome data segment
    :param label: label for the ratio barplot
    """
    plt.rcParams['mathtext.fontset'] = 'cm'
    # figfile = '/Users/davidmurphy/GoogleDrive/linked_selection/ls_paper/hipi-{}'.format(label.replace(' ', '-'))
    # replace 0 counts with nans
    hddhist[hddhist == 0] = np.nan
    gnhist[gnhist == 0] = np.nan
    plt.figure(figsize=(13, 7))
    plt.subplots_adjust(left=0.1, right=0.96, top=0.96)

    # plot of ratio
    ax = plt.subplot(3, 1, (2, 3))
    # normalize the data
    hddnorm = hddhist / np.sum(hddsites)
    gnnorm = gnhist / np.sum(gnsites)
    boundbox = dict(facecolor='White', alpha=0.8, edgecolor='black', boxstyle='round')
    ratio = r'$\mathrm{\frac{proportion\ of\ outliers}{proportion\ of\ num_bases}}$'
    plt.bar(left=np.arange(100), height=hddnorm / gnnorm, width=1, align='edge')
    plt.axhline(y=1.0, color='k', ls='--', lw=1, alpha=0.9)
    plt.text(0.5, 0.95, ratio, fontsize=24, ha='center', va='top', transform=ax.transAxes, bbox=boundbox)
    plt.ylabel('ratio', fontsize=16, labelpad=10)
    plt.yticks(fontsize=16)
    plt.xlabel(label, fontsize=16, labelpad=10)
    plt.xticks(fontsize=16)
    plt.xlim(-5, 105)

    # plot of neutral sites per bin
    plt.subplot(311)
    plt.plot(gnhist, 'o-', label='num_bases neutral sites')
    plt.plot(hddhist, 'o-', label='outlier neutral sites')
    plt.ylabel('neutral sites', fontsize=16, labelpad=10)
    plt.yticks(fontsize=16)
    plt.yscale('log')
    plt.xticks([0, 20, 40, 60, 80, 100], [])
    plt.xlim(-5, 105)
    plt.legend(prop={'size': 14})
    # plt.savefig(figfile)
    plt.show()

    return None


def plot_repeat_enrichment(enr_file):
    r = np.loadtxt(enr_file, dtype=str)[:, 1:].astype(float)
    plt.figure(figsize=(12, 7))
    message = []
    hipi_bases = np.sum(r[:, 1])
    all_bases = np.sum(r[:, 0])
    # plt.plot(np.cumsum(r[:, 0] / all_bases), 'o-', label='all')
    # plt.plot(np.cumsum(r[:, 1] / hipi_bases), 'o-', label='high')
    rclasses = 'Low_complexity Simple_repeat DNA LTR LINE SINE'.split()
    # a = 0
    for i in xrange(2, 13, 2):
        # plt.plot(np.cumsum(r[:, i]) / np.sum(r[:, i]), 'o-')
        # plt.plot(np.cumsum(r[:, i + 1]) / np.sum(r[:, i + 1]), 'o-')
        # plt.title(rclasses[a])
        # a += 1
        # plt.show()
        plt.subplot(211)
        plt.plot(range(1, 23), r[:, i + 1] / r[:, 1] * r[:, 0] / r[:, i], 'o-')
        plt.subplot(212)
        plt.plot(range(1, 23), np.log10(r[:, i + 1]), 'o-')
        # m = 'all: {:.4f}; high: {:.4f}'.format(np.sum(r[:, i]) / all_bases, np.sum(r[:, i + 1]) / hipi_bases)
        # message.append(m)

    # print '\n'.join('{}: {}'.format(c, m) for c, m in zip(rclasses, message))
    # print hipi_bases / all_bases
    plt.legend(rclasses)
    plt.xticks(range(1, 23))
    plt.xlim(0.5, 22.5)
    plt.grid(axis='both')
    plt.show()

    plt.figure(figsize=(12, 7))
    plt.bar(left=np.arange(6),
            height=np.sum(r[:, 3:14:2], axis=0) / hipi_bases * all_bases / np.sum(r[:, 2:13:2], axis=0))
    # plt.axhline(y=1, ls='--', color='Gray')
    plt.xticks(range(6), [st.replace('_', ' ') for st in rclasses])
    plt.grid(axis='y')
    plt.show()

    plt.figure(figsize=(12, 7))
    plt.plot(range(1, 23), np.cumsum(r[:, 1] / hipi_bases), 'o-', label='high')
    plt.plot(range(1, 23), np.cumsum(r[:, 0] / all_bases), 'o-', label='all')
    plt.xticks(range(1, 23))
    plt.yticks(np.arange(0, 1.01, 0.05))
    plt.grid(axis='both')
    plt.legend()
    plt.show()


def analyze_segments(npfile, init_file, repeat_dir):
    """
    A number of stats and plots for the segments containing high diversity data appearing in the tail of the
    sort-by-prediction plots. Different features of these segments are analyzed and compared to all segments of
    neutral data used in the inference. For example, the number of "callable" bases per segment using the 1000
    genomes mask file in high-pi vs. all segments, etc.
    :param npfile: file containing segments for high diversity data -- incl. all possible mapstruct arrays
    :param init_file: runstruct/mapstruct init file
    :param repeat_dir: a file containing repeat coordinates
    """
    vs = np.load(npfile)
    chroms = set(vs[:, 0])
    # rclasses = 'Low_complexity Simple_repeat DNA LTR LINE SINE'.split()
    # print '#chrom\ttotal_bases\thipi_bases\t' + '\t'.join('{s}-all\t{s}-hipi'.format(s=s) for s in rclasses)

    for ch in chroms:
        # chrom string
        chrom = 'chr{}'.format(int(ch))
        # chrom index
        ci = (int(ch) - 1,)
        # vs organization:
        # 0=chrom 1=start 2=end 3=n_sites 4=phast_score_sum 5=phast_bases 6=gc_count 7=ref_bases 8=hm_bases
        # 9=hm_subs 10=hom 11=het
        # data for high diversity segments in chromosome 'ch'
        hipi = vs[vs[:, 0] == ch]
        # sort by range start
        hipi = hipi[hipi[:, 1].argsort()]
        # get nsite counts
        hipi_nsites = hipi[:, 3]
        # segment ranges in high diversity
        hipi_ranges = hipi[:, 1:3].astype(int)
        # segment lengths
        hipi_segment_lengths = hipi_ranges[:, 1] - hipi_ranges[:, 0]

        # data for whole chrom
        mst = MapStruct(init=init_file, ci=ci, seg=True, div=True, con=True, gc=True)
        # mst.__init__(rst=mst)
        # return select arrays
        chrom_ranges, phast_data, chrom_nsites = mst.masked('ranges con nsites'.split())
        # mask for removing sites that are ALREADY in high pi sites for fishers test
        # hpmask = np.in1d(chrom_ranges[:, 0], hipi_ranges[:, 0], invert=True)
        # chrom_ranges, phast_data, chrom_nsites = [x[hpmask] for x in chrom_ranges, phast_data, chrom_nsites]
        # segment lengths
        chrom_segment_lengths = chrom_ranges[:, 1] - chrom_ranges[:, 0]
        # fraction of num_bases sites in repeats
        chrom_segment_sum = np.sum(chrom_segment_lengths)

        # callability fractions
        # pass_count(chrom_ranges, chrom_nsites, hipi_ranges, hipi_nsites, mst.call_mask[ci[0]])

        # aligned phastcons sites per segment
        # phast_sites = phast_data[:, 1]
        # # aligned site density per segment
        # phast_density = phast_sites / chrom_segment_lengths
        # # mean aligned site density
        # mean_phast_density = np.sum(phast_sites) / chrom_segment_sum

        hipi_segment_sum = np.sum(hipi_segment_lengths)
        # aligned phastcons site density
        # hipi_phast_density = hipi[:, 5] / hipi_segment_lengths
        # mean_hipi_phast_density = np.sum(hipi[:, 5]) / np.sum(hipi_segment_lengths)

        # # start the list of printouts
        # printout = [chrom, chrom_segment_sum, hipi_segment_sum]
        # # increase the starts by 1 so that segments don't overlap as required by intersect function
        # chrom_ranges[:, 0] += 1
        # # increase the starts by 1 so that segments don't overlap as required by intersect function
        # hipi_ranges[:, 0] += 1
        # for cls in rclasses:
        #
        #     # all repeats for this chromosome merged
        #     repeat_file = '{dr}/{cl}/{ch}.{cl}.hg19.bed.gz'.format(dr=repeat_dir, cl=cls, ch=chrom)
        #     with zopen(repeat_file, 'r') as f:
        #         repeats = [map(int, line.split('\t')[1:3]) for line in f if not line.startswith('#')]
        #     repeats = list(merge(repeats))
        #     # intersection of repeat segments and all data segments for the chrom
        #     chrom_repeat_intersect = np.array(list(intersect(repeats, chrom_ranges)))
        #     chrom_repeat_segments = chrom_repeat_intersect[:, 1] - chrom_repeat_intersect[:, 0]
        #     chrom_repeats_sum = np.sum(chrom_repeat_segments)
        #     # chrom_fraction_repeats = chrom_repeats_sum / chrom_segment_sum
        #     printout.append(int(chrom_repeats_sum))
        #
        #     # intersection with high diversity and repeat segments
        #     hipi_repeat_intersect = np.array(list(intersect(repeats, hipi_ranges)))
        #     hipi_repeat_segments = hipi_repeat_intersect[:, 1] - hipi_repeat_intersect[:, 0]
        #     # fraction of num_bases sites in repeats
        #     hipi_repeats_sum = np.sum(hipi_repeat_segments)
        #     # hipi_fraction_repeats = hipi_repeats_sum / hipi_segment_sum
        #     printout.append(int(hipi_repeats_sum))
        #
        # print '\t'.join(str(s) for s in printout)

        # # plots
        # plt.figure(figsize=(12, 7))
        # plt.hist(hipi_phast_density, bins=250, weights=hipi[:, 3], normed=1, label='fraction aligned', alpha=0.75)
        # plt.hist(phast_density, bins=250, weights=chrom_nsites, normed=1, label='fraction aligned (all)', alpha=0.75)
        #
        # # plot labels
        # t1 = '{}, n={:.2e}, mean={:.3f}'.format(chrom, np.sum(hipi[:, 5]), mean_hipi_phast_density)
        # t2 = '{}, n={:.2e}, mean={:.3f}'.format(chrom, np.sum(phast_sites), mean_phast_density)
        # plt.title(t1 + '\n' + t2)
        # plt.xlim(-0.1, 1.1)
        # # plt.title('{}, n={:.2e}'.format(chrom, np.sum(ar[:, 3])))
        # # plt.xlim(0, chromosome_length(chrom=chrom))
        # plt.legend()
        # plt.show()


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy'):
        # pass
        main_local()
    else:
        main_remote()
