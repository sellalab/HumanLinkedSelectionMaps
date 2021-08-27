import numpy as np

__author__ = 'davidmurphy'


def moving_avg(datafile, window=2.5e5, overlap=1.25e5, normalize=False):
    """
    A function that computes a moving average over data in a specified window (usually referring to spatial
    coordinate system of some kind) and a given overlap in the same coordinates.
    :param datafile: an array file of data that contains coordinates and associated values (.npy or .txt is OK)
    :param window: a window in the coordinate system of the data array (i.e. 1MB of data in a chromosome)
    :param overlap: the overlap between windows
    :param normalize: optional flag to normalize data values to 1
    :return bin_arr: a numpy.array of the format [window_start_position, avg_window_value]
    """

    if isinstance(datafile, str):

        # test if datafile is the name of a file or not
        if datafile.endswith('.npy'):
            data = np.load(datafile)
        else:
            data = np.loadtxt(datafile)
    else:
        # if not the name of a file, presumably it is an actual array
        data = datafile

    if normalize:
        data[:, 1] /= np.mean(data[:, 1])  # normalize data values

    # settings for numpy.searchsorted:
    # initial = lowest search bound, which yields index 0
    # final = highest search bound - overlap so that it will yield an index one step below the max position
    initial = data[:, 0].min()
    final = data[:, 0].max() - overlap
    increment = overlap
    search_points = np.arange(start=initial, stop=final, step=increment)

    # numpy.searchsorted(arr, x) will find the index of the value in arr that is the closest to x as long as the
    # value in arr is <= x (x can also be an array of values).
    # EXAMPLE:
    # numpy.searchsorted([0,1,2], [0.5, 1, 1.5])
    # array([1, 1, 2])
    low_idx = np.searchsorted(data[:, 0], search_points)
    high_idx = np.searchsorted(data[:, 0], search_points + window)

    assert len(low_idx) == len(high_idx)

    rows = len(low_idx)
    bin_arr = []
    # use a cutoff to throw away any windows that have less than 5% of the average count of data points
    cutoff = np.mean(high_idx - low_idx) * 0.05

    for i in xrange(rows):

        start = low_idx[i]  # get indices for start and end of each row
        end = high_idx[i]
        pos = data[start, 0]  # use the lower bound of each window for position, i.e. data[start]

        if end - start > cutoff:

            # skip windows with very little data below the cutoff count
            # take the mean of the values in each window of data if count is above cutoff
            # add to nested array-like list bin_arr
            bin_avg = np.mean(data[start:end, 1])
            bin_arr += [[pos, bin_avg]]

        else:
            # record None for sites below cutoff (better for plotting)
            bin_arr += [[pos, None]]

    bin_arr = np.array(bin_arr)  # convert to array for easy indexing
    return bin_arr

# DEBUG:
# import project_dirs as pd
# ch = 1
# # pi  = '{d}compare/chr{c}_neut_pi.npy'.format(d=pd.stored_results, c=ch)
# div = '{d}chr{c}_h2m_div_20kb_windows.txt'.format(d=pd.mutprox, c=ch)
# # a25fm = '{d}compare/AA_Map_downSample_25pct_cmmb_0_neutral_mutrate_chr1.LS'.format(d=pd.stored_results, c=ch)
# a15m = '{d}compare/AA_Map_downSample_15pct_neutral_mutrate_chr1.LS'.format(d=pd.stored_results, c=ch)
# a15m_7t = '{d}compare/AA_Map_downSample_15pct_neutral_mutrate_chr1_7t.LS'.format(d=pd.stored_results, c=ch)
# #
# #
# bmap7 = np.loadtxt(a15m_7t)
# bmap3 = np.loadtxt(a15m)
#
# mut = np.loadtxt(div)
#
# correction = np.interp(bmap7[:,0], mut[:,0], mut[:,1])
#
# bmap7[:,1] *= correction  # multiply by divergence
# bmap3[:,1] *= correction
# #
# # # a1_7m = '{d}compare/AA_Map_downSample_1pct_neutral_mutrate_chr1.LS'.format(d=pd.stored_results, c=ch)
# # # aa_m  = '{d}trial 63/AA_Map_downSample_25pct_cmmb_0_neutral_mutrate_chr{c}.LS'.format(d=pd.stored_results, c=ch)
# # # aa_c = '{d}trial 64/AA_Map_downSample_25pct_cmmb_0_neutral_constant_chr{c}.LS'.format(d=pd.stored_results, c=ch)
# # # dc_m = '{d}trial 65/deCODE_downSample_25pct_cmmb_0_neutral_mutrate_chr{c}.LS'.format(d=pd.stored_results, c=ch)
# # # dc_c = '{d}trial 66/deCODE_downSample_25pct_cmmb_0_neutral_constant_chr{c}.LS'.format(d=pd.stored_results, c=ch)
# # # factor = np.interp(a15[:,0], div[:,0], div[:,1])
# # # d25fm = '{d}compare/deCODE_downSample_25pct_cmmb_0_neutral_mutrate_chr1.LS'.format(d=pd.stored_results, c=ch)
# # # d15m = '{d}compare/deCODE_downSample_15pct_neutral_mutrate_chr1.LS'.format(d=pd.stored_results, c=ch)
# # # d25f = '{d}compare/deCODE_downSample_25pct_cmmb_0_neutral_constant_chr1.LS'.format(d=pd.stored_results, c=ch)
# # # d15 = '{d}compare/deCODE_downSample_15pct_neutral_constant_chr1.LS'.format(d=pd.stored_results, c=ch)
# #
# # a = mv_avg(pi, 2.5e5, 1.25e5)
# b = moving_avg(bmap3, 2.5e5, 1.25e5)
# e = moving_avg(bmap7, 2.5e5, 1.25e5)
# # f = moving_avg(a1_7m, 2.5e5, 1.25e5)
# # #
# # # d =np.array(mv_avg(div, 2.5e5, 1.25e5), dtype='f8')
# # # v = mv_avg(d25f, 2.5e5, 1.25e5)
# # # w = mv_avg(d15, 2.5e5, 1.25e5)
# # # x = mv_avg(d25fm, 2.5e5, 1.25e5)
# # # y = mv_avg(d15m, 2.5e5, 1.25e5)
# #
# # plt.title('7t plot')
# # # plt.plot(factor[:,0], factor[:,1])
# # plt.plot(a[:,0], a[:,1], alpha=0.5, linewidth=2, color='k')
# plt.plot(b[:,0], b[:,1], alpha=0.5, linewidth=2, color='c')
# # # plt.plot(c[:,0], c[:,1], alpha=0.5)
# # # plt.plot(d[:,0], d[:,1], alpha=0.5)
# plt.plot(e[:,0], e[:,1], alpha=0.5, linewidth=2, color='r')
# # # plt.plot(f[:,0], f[:,1], alpha=0.5)
# #
# # # plt.plot(v[:,0], v[:,1], alpha=0.5)
# # # plt.plot(w[:,0], w[:,1], alpha=0.5)
# # # plt.plot(x[:,0], x[:,1], alpha=0.5)
# # # plt.plot(y[:,0], y[:,1], alpha=0.5)
# #
# # plt.legend(["pi", "b3", "b7"])
# #
# plt.show()
