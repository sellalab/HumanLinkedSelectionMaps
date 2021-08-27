from sys import argv

from cluster_code.inference.results.mapfunctions import *

__author__ = 'davidmurphy'


def main_local():
    import matplotlib.pyplot as plt
    plt.figure(figsize=(16, 9))
    leg = ['prPi', 'prPred', 'mcPred']
    plotpi = False
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/'
    for f in 'primateCons95segs_BS1_CS0_161108110259_final', 'mcvicker_map_BS1_CS0_161129211053_initial':
        f = rdir + f + '.txt'
        rst = RunStruct(init=f)
        ch = 'chr1'
        scale = 1e6
        pos, pi, pred = map_diversity(rst=rst, chrom=ch, scale=scale).T
        if not plotpi:
            plt.plot(pos, pi, lw=3, color='k', alpha=0.75)
            plotpi = True
        plt.plot(pos, pred, lw=2, alpha=0.5)
    plt.legend(leg)
    plt.show()


def main_remote():
    if len(argv) == 5:
        init_file = argv[1]
        chrom = argv[2]
        scale = float(argv[3])
        out_file = argv[4]
        rst = RunStruct(init=init_file)
        map_array = map_diversity(rst=rst, chrom=chrom, scale=scale)
        np.save(file=out_file, arr=map_array)
    else:
        print 'usage: map_diversity <init_file> <chrom> <scale> <out_file>'


def map_diversity(rst, chrom, scale):
    """
    This function is used to calculate rsquared for observed vs. predicted diversity levels at a given spatial scale
    and for a given set of inferred parameters and data. The scale-rsquared values are printed to stderr and also
    sent to a file.
    :param rst: a RunStruct containing the params
    :param chrom: human chromosome
    :param scale: scale size in bp on which to calculate rsquared
    :return [pos, pi, pred] array
    """
    # data vectors (positions, hom, het, u, predictions)
    vectors = infmap_chrom(rst=rst, chrom=chrom)
    # get the mean predicted and observed diversity levels for the current spatial scale
    bins, pi, predicted = binned_map(vectors=vectors, binsize=scale, stepsize=0.5 * scale)

    return np.column_stack([bins, pi, predicted])


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
    else:
        main_remote()
