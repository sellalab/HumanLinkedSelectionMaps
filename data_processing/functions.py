import os
import re
import numpy as np
from scipy.stats import norm
from datetime import timedelta
from gzip import open as zopen

__author__ = 'davidmurphy'


"""
General functions for reference tables, basic statistics and string formatting
"""


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
       
    WEB REFERENCE: http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    """
    # import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def harmonic_mean(values):
    """return the harmonic mean of a set of values"""
    den = len(values)
    num = 0
    for v in values:
        num += 1.0 / v

    return (num / den) ** -1.0


def rsquared_function(xobs, yobs):
    """
    Calculates the coefficient of determination from two vectors of observations and predictions.
    Adapted from: https://en.wikipedia.org/wiki/Coefficient_of_determination
    :param xobs: the data that should be predicted by y = x
    :param yobs: the predictions at each point for x
    :return cd: the coefficient of determination
    """
    # only pass cleaned data to this function
    assert (np.isfinite(xobs) & np.isfinite(yobs)).all()
    assert len(xobs) == len(yobs)
    x_bar = np.mean(xobs)  # mean value of observed data
    ss_tot = np.sum(np.square(xobs - x_bar))  # total sum of squares -- proportional to variance in xobs data
    ss_res = np.sum(np.square(xobs - yobs))  # total sum of residuals
    rsq = 1 - ss_res / ss_tot  # calculate R^2
    return rsq


def rsquared_residuals(xobs, yobs):
    # only pass cleaned data to this function
    assert (np.isfinite(xobs) & np.isfinite(yobs)).all()
    assert len(xobs) == len(yobs)
    x_bar = np.mean(xobs)  # mean value of observed data
    y_bar = np.mean(yobs)
    ss_tot = np.sum(np.square(xobs - x_bar))  # total sum of squares -- proportional to variance in xobs data
    residual = 2.0 * np.sum(np.square(xobs - yobs) * (yobs - y_bar)) / ss_tot

    return residual


def errorprop_fraction(mx, my, vx, vy):
    """
    Calculate propogation of error for random variable Z, where Z = Y/X and X & Y are indepedent random variables with
    means (mx, my) and variances (vx, vy), respectively
    """
    return mx**-2.0 * ((vx * (my / mx)**2.0) + vy)


def relative_error(n, k):
    """
    Calculate the relative error for a sample of data assumed to follow a binomial distribution
    :param n: number of trials
    :param k: successful trials
    :return relative_err: the relative error of the sample
    """
    # convert n & k from integers to floats or arrays of floats
    if isinstance(n, int) or isinstance(k, float):
        n = float(n)
        k = float(k)
    else:
        assert isinstance(n, np.ndarray) and isinstance(k, np.ndarray)
        n = n.astype('f8')
        k = k.astype('f8')

    # rate estimate (p) for a random variable X ~ Binom(n, p): E[X/n] -> E[X]/n
    phat = k / n
    qhat = 1 - phat
    # standard error for rate estimate (p): sqrt{Var[X/n]}
    # -> sqrt{Var[X]/n^2}
    # -> sqrt{Var[X]}/n)
    stderr = np.sqrt(phat * qhat / n)
    # relative error is the standard error scaled by the rate estimate
    relative_err = stderr / phat

    return relative_err


def wilcoxon_signed(x1, x2):
    """
    Wilcoxon signed-rank test
    =========================
    Assumptions
    -----------
    + Data are paired and come from the same population.
    + Each pair is chosen randomly and independently.
    + The data are measured at least on an ordinal scale (i.e., they cannot be nominal).
    Test procedure
    --------------
    + Let n be the sample size, i.e., the number of pairs. Thus, there are a total of 2n data points.
    + For pairs i = 1,...,n, let x1,i and x2,i denote the measurements.
    + H0: difference between the pairs follows a symmetric distribution around zero
    + H1: difference between the pairs does not follow a symmetric distribution around zero.
        1. For i = 1,...,n, calculate abs(x2,i - x1,i) and sign(x2,i - x1,i)
        2. Exclude i and reduce n by 1 (n->nr) for each i where abs(x2,i - x1,i) == 0
        3. Order the remaining nr pairs from smallest absolute difference to largest absolute difference
        4. Rank the pairs, starting with the smallest as 1. Ties receive a rank equal to the average of
        the ranks they span. Let ri denote the rank.
        5. Calculate the test statistic w:
            w = sum(sign(x2,i - x1,i) * ri) for i in nr
        6. Under null hypothesis, w follows a specific distribution with no simple expression. This
        distribution has an expected value of 0 and a variance:
            var_w = (nr / 6) * (nr + 1) * (2 * nr + 1)
        w can be compared to a critical value from a reference table. The two-sided test consists in
        rejecting H0 if abs(w) >= w_critical_nr
    + As nr increases, the sampling distribution of w converges to a normal distribution.
    + Thus, for nr >= 10, a z-score can be calculated as z = (w / sigma_w); sigma_w = sqrt(var_w)
    + The two-sided test consists in rejecting H0 if z > z_critical_nr
    (source: https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test)
    :param x1: vector of numerical observations
    :param x2: vector of numerical observations
    :type x1: numpy.ndarray
    :type x2: numpy.ndarray
    :return p-value:
    """
    assert len(x1) == len(x2)
    # 1) abs difference, sign difference x1, x2
    dx = x2 - x1
    # sort dx by abs difference --> ranks = range(1, nr+1)
    dx = dx[abs(dx).argsort()]
    # get the signs *sorted* dx
    sx = np.sign(dx)
    # 2) mask zero difference indices
    if any(dx == 0):
        dx, sx = dx[(dx != 0)], sx[(dx != 0)]
    nr = len(dx)
    # 3) - 4) rank the differences from 1 to nr
    ri = np.arange(1, nr + 1)
    # deal with ties
    if len(np.unique(abs(dx))) < len(dx):
        pass
    # count sites where x2 > x1 and vice versa
    nx2 = np.sum(sx == 1, dtype='f8')
    nx1 = np.sum(sx == -1, dtype='f8')
    assert nx1 + nx2 == nr
    # 5) calculate test statistic
    # r_plus = np.sum(s_plus * ri)
    # r_minus = np.sum(s_plus * ri)
    # w = r_plus - r_minus
    w = np.sum(sx * ri)
    var_w = (nr / 6.0) * (nr + 1.0) * (2.0 * nr + 1.0)
    # 6) calculate p-value
    if nr < 10:
        print 'Warning! small sample size -- nr = {}'.format(nr)
    z = (w / np.sqrt(var_w))
    # use survival function to get 2-tailed p-value with abs(z)
    pval = 2.0 * norm.sf(abs(z))
    frac = nx2 / nr
    # return na / (na+nb), p_value
    return frac, pval


def count_all_pairs(sample_size, alt_count):
    """
    Mini-function to calculate the number of monomorphic and polymorphic pairs in a sample of sites with ref/alt counts
    :param sample_size: total number of chromosomes in the sample
    :param alt_count: total number of non-reference chromosomes in the sample
    :return mono_pairs, poly_pairs: a tuple of the counts for each type of pairing
    """
    # cast inputs as numpy.array so that numpy.sum function is sure to work on results
    sample_size = np.array(sample_size)
    alt_count = np.array(alt_count)
    total_pairs = 0.5 * sample_size * (sample_size - 1)
    poly_pairs = alt_count * (sample_size - alt_count)
    mono_pairs = total_pairs - poly_pairs

    return mono_pairs, poly_pairs


def calc_pi(sample, alts):
    # Calculate pi from sample, alternate counts
    return 2.0 * alts * (sample - alts) / (sample * (sample - 1))


def pairwise_divergence(r1, a1, r2, a2, altmatch=True):
    """
    Calculate the pairwise divergence between two population samples.
    r1, a1, r2, a2 are reference and alternate counts for populations
    1 and 2, respectively.
    If altmatch is False, then the allele at a1 does not equal the 
    allele at a2 and and additional calculation is done for fi
    :return avg_div: average pairwise divergence between two pops 
    """
    # convert everything to floating point for safety (works on arrays and single numbers)
    r1, a1, r2, a2 = [1.0 * x for x in r1, a1, r2, a2]

    n1 = r1 + a1  # population 1 sample
    n2 = r2 + a2  # population 2 sample
    nc = n1 * n2  # total possible pairwise comparisons between populations
    fi = (r1 * a2) + (a1 * r2)  # indicator function fi = [1 for mismatches, 0 for matches]

    # additional comparison for mismatching alternate alleles
    if not altmatch:
        fi += (a1 * a2)

    # the fraction of mismatching comparisons over possible comparisons gives the average pairwise divergence
    return fi / nc


def weighted_variance(x, w):
    """
    Returns the variance for weighted data x (e.g., weighted variance)
    :param x: array of data 
    :param w: array of weights
    """
    # first get weighted mean with weighted average
    xmean = np.average(x, weights=w)
    # use weighted average again to return the weighted average of the variance denominator: (x-xmean)^2
    return np.average(np.power(x - xmean, 2), weights=w)


def scientific_pow10(num):
    # return the leading float and power of 10 that goes with it
    sci_string = '{:.18e}'.format(num)
    lead_float, pow10 = sci_string.split('e')
    return float(lead_float), int(pow10)


def expocdf(mu, x):
    """
    returns the area between [0, x] for an exponential with parameter mu, i.e. the cumulative density function of the
    expo(mu).
    :param mu: expo rate parameter. for the pdf of expo(mu): mean = 1/mu, variance = 1/mu^2
    :param x: upper lim of the integral
    :return cdf(x): the area under expo(mu) from [0, x]
    """
    return 1 - np.exp(-x * mu)


def chromosome_length(chrom):
    # return the integer chromosome length in bp for chrom
    chr_length = dict(chrY=59373566, chrX=155270560, chr13=115169878, chr12=133851895, chr11=135006516, chr10=135534747,
                      chr17=81195210, chr16=90354753, chr15=102531392, chr14=107349540, chr19=59128983, chr18=78077248,
                      chr22=51304566, chr20=63025520, chr21=48129895, chr7=159138663, chr6=171115067, chr5=180915260,
                      chr4=191154276, chr3=198022430, chr2=243199373, chr1=249250621, chr9=141213431, chr8=146364022,
                      chr2L=23011544, chr2R=21146708, chr3L=24543557, chr3R=27905053)
    return chr_length[chrom]


def time2str(time_obj):
    # use a strict hours:minutes:seconds format to store time durations
    assert isinstance(time_obj, timedelta)
    tot_sec = int(time_obj.total_seconds())
    seconds = tot_sec % 60
    minutes = (tot_sec % 3600) / 60
    hours = tot_sec / 3600
    return '{}:{}:{}'.format(hours, minutes, seconds)


def str2time(time_str):
    # convert "n day(s), hours:minutes:seconds" to a datetime.timedelta object
    assert isinstance(time_str, str)
    retime = re.compile(r'((?P<days>\d+)\sdays*,\s)*(?P<hours>\d{1,2}):(?P<minutes>\d{1,2}):(?P<seconds>[\d.]+)')
    try:
        timedict = retime.search(time_str).groupdict()
        for (k, v) in timedict.items():
            if v:
                timedict[k] = float(v)
            else:
                del timedict[k]
        return timedelta(**timedict)
    except ValueError:
        raise ValueError('Elapsed time format <{}> not recognized'.format(time_str))
    except AttributeError:
        raise AttributeError('Elapsed time format <{}> not recognized'.format(time_str))


def swap_root(filepath, other_strings=None):
    """
    A python implementation of "swap_root.awk". This function opens a file and changes the root directory to the
    current environment if it does not already match. The file is saved under the same name.
    :param filepath: an init file path
    :param other_strings: optional list of (string1, string2) tuples to swap
    :type filepath: str
    """
    # root paths for my computer and my cluster account
    local_root = '/Users/davidmurphy/GoogleDrive/linked_selection'
    clust_root = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection'
    # read in the content of the file as one string
    with open(filepath, 'r') as f:
        content = f.read()

    # swap the paths throughout the file
    if os.getcwd().startswith('/Users/davidmurphy'):
        new_content = content.replace(clust_root, local_root)
    else:
        new_content = content.replace(local_root, clust_root)

    # swap additional strings if specified
    if other_strings is not None:
        for (string1, string2) in other_strings:
            new_content = new_content.replace(string1, string2)

    # re-write file with updated root directory for current environment
    with open(filepath, 'w') as f:
        f.write(new_content)
