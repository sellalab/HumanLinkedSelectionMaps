from sys import argv
from classes.mapstruct import MapStruct, np, os
from classes.runstruct import RunStruct
from scipy.stats import pearsonr as pcor
import matplotlib.pyplot as plt
import seaborn

__author__ = 'davidmurphy'


def main_remote():
    if len(argv) == 3:
        init_file = argv[1]
        num_bins = eval(argv[2])
        rst = RunStruct(init=init_file)
        arr_file = '{}/result/sort_map/{}.{}.sorted.{}bins.npy'
        arr_file = arr_file.format(rst.root, rst.neut, rst.label, num_bins)
        arr = basic_sort(init=init_file, num=num_bins)
        np.save(arr_file, arr)


def main_local():
    num_bins = 100
    rdir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/final_files/'
    # folder = 'alternateclean'
    folder = 'phyloclean'
    # folder = 'pr90to99clean'
    d = rdir + folder

    files = ['{}/{}'.format(d, f) for f in os.listdir(d)]

    # titles = ['{}%'.format(p) for p in xrange(90, 100)]
    # titles = 'CADD McVicker Genic/Nongenic-cons'.split()
    titles = [f.split('/')[-1].split('.')[0] for f in files]

    for (f, t) in zip(files, titles):
        plot_result(init_file=f, folder=folder, title=t)
    plt.show()


def plot_result(init_file, folder, title, num_bins=100):
    from functions import swap_root
    swap_root(init_file)
    rst = RunStruct(init=init_file)
    arr_file = '{}/result/sort_map/{}/{}.{}.sorted.{}bins.npy'.format(rst.root, folder, rst.neut, rst.label, num_bins)
    figdir = '/Users/davidmurphy/GoogleDrive/linked_selection/murphy_et_al/figures/{}'.format(folder)
    if not os.path.isdir(figdir):
        os.mkdir(figdir)
    fig_file = '{}/{}.png'.format(figdir, rst.label)

    # rst.label = '{}.BS{}.{}.CS{}.{}'.format(rst.token, rst.bnum, rst.bsgrid, rst.cnum, rst.csgrid)
    # arr_file = '{}/result/sort_map/primate90to99/{}.{}.sorted.{}bins.npy'.format(rst.root, rst.neut, rst.label,
    #                                                                              num_bins)
    # arr = basic_sort(rst=rst, num=num_bins)
    # np.save(arr_file, arr)
    arr = np.load(arr_file)

    # dv, pi, pr, ar1, ar2, hnd = arr.T
    dv, pi, pr = arr.T

    meanpi = rst.stat.meanpi / rst.stat.meanpi0
    x = np.arange(len(dv))
    pi /= rst.stat.meanpi0
    pr /= rst.stat.meanpi0
    # ar1 /= ar1.mean()
    # ar1 *= meanpi
    # ar2 /= ar2.mean()
    # ar2 *= meanpi
    # hnd /= hnd.mean()
    # hnd *= meanpi

    corr, pval = pcor(pi, pr)

    gray = 'DarkSlateGray'
    fig = plt.figure(figsize=(10, 7))
    plt.subplots_adjust(left=0.13, bottom=0.12, right=0.78, top=0.90, wspace=0.05, hspace=0.2)
    fig.add_subplot(111)

    plt.plot(x, pi, label='observed diversity', color=gray, lw=3)
    plt.plot(x, pr, label='predicted', color='Fuchsia', alpha=0.75, lw=3)

    # plt.plot(x, ar1, label='archaic haplotype probability (Song)', color='DodgerBlue', alpha=0.75, lw=3)
    # plt.plot(x, ar2, label='archaic haplotype probability (Reich)', color='Purple')
    # plt.plot(x, hnd, label='scaled YRI-Altai divergence', color='FireBrick', alpha=0.75, lw=3)

    message = ' {}\n {:>13}'.format(*'no background;selection'.split(';'))
    plt.axhline(y=1, color='k', ls='--', lw=2.0)
    # plt.text((0.5 * num_bins), 1.005, 'w/o background selection', ha='center', va='bottom', fontsize=18)
    plt.text((1.025 * num_bins), 1, message, ha='left', va='center', fontsize=20)

    message = '   {:>6}\n   {}'.format(*'mean;diversity'.split(';'))
    plt.axhline(y=meanpi, color=gray, ls=':', lw=2.0)
    # plt.text((0.8 * num_bins), meanpi - 0.005, 'mean diversity', ha='center', va='top', fontsize=18, color=gray)
    plt.text((1.025 * num_bins), meanpi, message, ha='left', va='center', fontsize=20, color=gray)

    # Pearson correlation
    plt.text((0.5 * num_bins), 0.55, r'$Pearson\ R^2 = {:.4f}$'.format(corr), ha='center', va='center', fontsize=20)

    # title w/ species or cons info
    plt.title(title, fontsize=20)
    # plt.title(rst.token, fontsize=20)
    # plt.title(init_file.split('/')[-1].split('.')[0], fontsize=20)

    plt.xlabel('strength of background selection', fontsize=20)
    plt.xticks(fontsize=18)
    plt.xlim((-0.01 * num_bins), (1.01 * num_bins))

    plt.ylabel('scaled diversity', fontsize=20, labelpad=10)
    plt.yticks(fontsize=18)
    plt.ylim(0.37, 1.23)

    plt.legend(prop={'size': 18}, loc='upper left', ncol=1, numpoints=3, borderaxespad=0.1)
    plt.savefig(fig_file)
    # exit(1)

    # arr = basic_sort(rst=rst, num=num_bins)
    # np.save(arr_file, arr)


def basic_sort(init, num):
    """
    Sort data arrays by predictions and average results into bins.
    :type init: RunStruct innit file
    :param num: the number of bins to sort into
    """
    # load basic arrays and get prediction sory index
    mst = MapStruct(init=init, seg=True, div=True)
    sidx = mst.prediction_sorting_index

    # sort predictions and basic div/poly data
    pred = mst.prediction(umap=mst.uconst)[sidx].astype('f8')
    nt, div = [a[sidx].astype('f8') for a in mst.nt, mst.div]
    sites, subs = div.T

    # get start and end indices per partition
    idx = sortbin_edges(sites=sites, numbins=num)

    # gather data summaries
    sorted_array = []
    for (i, j) in idx:

        # calculate each statistic for the present bin
        norm_div = np.sum(subs[i:j]) / (mst.stat.meandiv * np.sum(sites[i:j]))
        scaled_pi = np.sum(nt[i:j, 1]) / (np.sum(nt[i:j]) * norm_div)
        # use weighted mean for pred
        prediction = np.average(pred[i:j], weights=sites[i:j])
        sorted_array.append([norm_div, scaled_pi, prediction])

    return np.array(sorted_array)


def sortbin_edges(sites, numbins):
    """get the upper indices of sorted data array that divide data into bins of equal neutral site counts"""

    # get the number of sites needed in each bin such that numbins x numsites = total sites
    numsites = int(np.sum(sites) / numbins)

    # find indices that partition cumulative sorted site count into numsites
    cumsites = np.cumsum(sites)
    bounds = np.arange(numsites, cumsites[-1], numsites)

    # get the ending index for each partition
    jdx = list(np.searchsorted(a=cumsites, v=bounds))

    # return a list of (start, end) indices for each partition
    return zip([0] + jdx[:-1], jdx)


if __name__ == '__main__':
    if os.getcwd().startswith('/Users/davidmurphy/'):
        main_local()
        # other_sorts(500)
    else:
        main_remote()
