__author__ = 'davidmurphy'

import os
import gzip
import seaborn
import numpy as np
import matplotlib as mpl
# from skmisc import loess
from sys import stderr, stdout
from neutral_data import FilePaths, get_neutral_masks, token
from classes.runstruct import ChromStruct
from data_processing.data_tools import randint_unique

# set high DPI on figures globally
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 12

plt = mpl.pyplot


"""
BINNED DATA ANALYSIS TOOLS
"""


class BinnedData:
    """a structure that holds binned data and its properties"""
    def __init__(self, ch, neut, outg, sortby):
        # initialize FilePaths object to retrieve saved files
        fpth = FilePaths()

        # load bins edges and bin counts
        za = np.load(fpth.fbincnt(sortby))
        self.bins = za['bins']
        self.cnt = za[ch+'cnt']  # NOTE: ch is either x or a for X/Aut

        # load pi & div
        self.pi = np.load(fpth.fbinpi(neut, sortby))[ch+'pi']
        self.div = np.load(fpth.fbindiv(outg, sortby))[ch+'div']

        # create an x-axis label for plots baesd on the sorting type
        if sortby == 'dist':
            self.xlab = 'cM to nearest exon'
        else:
            self.xlab = 'B values'

    @property
    def msk(self):
        """mask to remove sites where pi=nan or cnt=0 or div=nan or div=0"""
        return np.isfinite(self.pi * self.div) & (self.cnt * self.div != 0)

    @property
    def mskbins(self):
        """masked bins"""
        return self.bins[self.msk]

    @property
    def mskcnt(self):
        """masked counts"""
        return self.cnt[self.msk]

    @property
    def mskpi(self):
        """masked pi values"""
        return self.pi[self.msk]

    @property
    def mskdiv(self):
        """masked D values"""
        return self.div[self.msk]

    @property
    def pi_div(self):
        """pi/D (use masked pi & div by default to avoid division by 0)"""
        return self.mskpi / self.mskdiv


def fit_loess(xi, yi, wts, span):
    """get loess fitted values for each xi input"""
    from skmisc import loess
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()

    return lo.outputs.fitted_values


def predict_loess(xi, yi, wts, span, xtest):
    """get loess predicted values at new points from xtest"""
    from skmisc import loess
    lo = loess.loess(xi, yi, weights=wts, span=span)
    lo.fit()
    predict = lo.predict(xtest)

    # NOTE: copy prevents a bug where repeated calls affected previous results
    return np.copy(predict.values)


def plot_1(a, x, span):
    """plot pi vs cM to exon or B value"""
    # plot log10 neutral site counts per bin
    plt.figure(figsize=(10, 6))
    plt.subplot(311)
    plt.scatter(a.bins, np.log10(a.cnt), color='mediumorchid')
    plt.scatter(x.bins, np.log10(x.cnt), color='darkturquoise')
    plt.ylabel('log10 sites')
    plt.xticks(color='white')

    # establish shared x limits 2.5% beyond width of bin range
    offset = (a.bins.max() - a.bins.min()) * 0.025
    xmin = a.bins.min()-offset
    xmax = a.bins.max()+offset
    plt.xlim(xmin, xmax)

    # plot LOESS smoothed pi for aut
    plt.subplot(3, 1, (2, 3))
    lbl = r'$\bar{\pi_A}$'
    clr = 'mediumorchid'
    api = predict_loess(a.mskbins, a.mskpi, a.mskcnt, span, a.bins)
    plt.plot(a.bins, api, color=clr, label=lbl)

    # plot LOESS smoothed pi for X
    xpi = predict_loess(x.mskbins, x.mskpi, x.mskcnt, span, x.bins)
    lbl = r'$\bar{\pi_X}$'
    clr = 'darkturquoise'
    plt.plot(x.bins, xpi, color=clr, label=lbl)

    # label and show the plot
    plt.xlabel(a.xlab)
    plt.xlim(xmin, xmax)
    plt.ylabel(r'$\bar{\pi}$')
    plt.legend(loc='upper center')
    plt.show()


def plot_2(a, x, span):
    """plot div vs cM to exon or B value"""
    # plot LOESS smoothed div for aut
    plt.figure(figsize=(10, 6))
    lbl = r'$\bar{D_A}$'
    clr = 'mediumorchid'
    adiv = predict_loess(a.mskbins, a.mskdiv, a.mskcnt, span, a.bins)
    plt.plot(a.bins, adiv, color=clr, label=lbl)

    # plot LOESS smoothed div for X
    xdiv = predict_loess(x.mskbins, x.mskdiv, x.mskcnt, span, x.bins)
    lbl = r'$\bar{D_X}$'
    clr = 'darkturquoise'
    plt.plot(x.bins, xdiv, color=clr, label=lbl)

    plt.xlabel(a.xlab)
    plt.ylabel(r'$\bar{D}$')
    plt.legend()
    plt.show()


def plot_3(a, x, span):
    """plot pi/div vs cM to exon or B value"""
    # plot LOESS smoothed pi/div for aut
    plt.figure(figsize=(10, 6))
    apd = predict_loess(a.mskbins, a.pi_div, a.mskcnt, span, a.bins)
    lbl = r'$\bar{\pi_A}/\bar{D_A}$'
    clr = 'mediumorchid'
    plt.plot(a.bins, apd, color=clr, label=lbl)

    # plot LOESS smoothed div for X
    xpd = predict_loess(x.mskbins, x.pi_div, x.mskcnt, span, x.bins)
    lbl = r'$\bar{\pi_X}/\bar{D_X}$'
    clr = 'darkturquoise'
    plt.plot(x.bins, xpd, color=clr, label=lbl)

    plt.xlabel(a.xlab)
    plt.ylabel(r'$\bar{\pi}/\bar{D}$')
    plt.legend()
    plt.show()


def plot_4(a, x, span):
    """plot X:A pi/div vs cM to exon or B value"""
    # get pi/div for X/A
    apd = predict_loess(a.mskbins, a.pi_div, a.mskcnt, span, a.bins)
    xpd = predict_loess(x.mskbins, x.pi_div, x.mskcnt, span, x.bins)

    plt.figure(figsize=(10, 6))
    plt.plot(a.bins, xpd/apd, color='forestgreen')

    plt.xlabel(a.xlab)
    lbl = r'$\frac{\bar{\pi_X}/\bar{D_X}}{\bar{\pi_A}/\bar{D_A}}$'
    plt.ylabel(lbl, fontsize=16)
    plt.show()


def main1():
    neut = 'YRI'
    outg = 'rheMac3'
    sortby = 'bval'
    span = 0.1
    atab = BinnedData('a', neut, outg, sortby)
    xtab = BinnedData('x', neut, outg, sortby)

    plot_1(atab, xtab, span)
    plot_2(atab, xtab, span)
    plot_3(atab, xtab, span)
    plot_4(atab, xtab, span)
