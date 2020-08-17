# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import pandas as pd
from verkko.binner import bins as binner
from scipy.stats import binned_statistic, binned_statistic_2d, rankdata
from scipy.optimize import leastsq
from scipy.ndimage import convolve

def plot_powerlaw(values, main='', xlabel='', ylabel='', fig=None, ax=None, label=''):
    """
    Function for plotting power-log distribution. If many plots are to be added, iterate with fig and ax.
    (values): list or 1d array of values to plot distribution
    (main): plot title
    (xlabel): value
    (ylabel):
    (fig): pyplot figure if plotting various graphs
    (ax): pyplot subplot if plotting various graphs
    (label): name of the values to display
    """
    if not fig:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(main)
    n_val = len(values)
    if not isinstance(values[0], int):
        values = [int(v) for v in values]
    p = np.bincount(values)/float(n_val)
    x = np.linspace(0, max(values), len(p))
    ax.loglog(x, p, '.', label=label)
    ax.legend(loc=0)
    return fig, ax

def powerlaw(values, main='', xlabel='', ylabel='', label='', fig=None, ax=None, factor=1.2, fit=True):

    """
Plot powerlaw distribution by binning values into log-bins defined by factor. This version also fits a powerlaw distribution of the form:
    P(x) = a*exp(-x/xc)(x+x0)**(-alpha)

Parameters:
    (values): list or 1D array of values to obtain a density from
    (main): matplotlib title
    (xlabel): x-axis label
    (ylabel): y-axis label
    (label): label to appear on the graph
    (fig): matplotlib figure
    (ax): matplotlib axis
    (factor): factor > 1 to log-bin values using verkko.

Returns:
    (fig, ax): matplotlib figure and axis
    (p1): vector with fit coefficientes: [log(a), alpha, x0, xc]

    """
    if not fig:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(main)

    n_val = len(values)
    start = np.min(values)
    end = np.max(values)
    bins = binner.Bins(float, int(min(values)) + 1, max(values) + 1, 'log', factor)

    ind = bins.widths > 1
    counts, _,  _ = binned_statistic(values, values, statistic='count', bins=bins.bin_limits)
    counts = counts[ind]/float(n_val)
    bins_c = bins.centers[ind]
    ax.loglog(bins_c, counts, '.',  label=label)

    try:
        if fit:
            p1, succ = powerlaw_coeff(bins_c, counts)
            y_fit = 10**(powerlaw_func(p1, bins_c))
            ax.loglog(bins_c, y_fit, label= 'Fit: ' + r'$\gamma = $' + str(round(p1[1], 2)), color='blue')
        else:
            p1 = None
    except:
        p1 = None

    ax.legend(loc=0)
    return fig, ax, p1


powerlaw_func = lambda p, x: p[0] + -p[1]*np.log(x + p[2]) -x/p[3]
errfunc = lambda p, x, y: (powerlaw_func(p, x) - y)**2

def powerlaw_coeff(x, y):

    p0 = [0.9, 1.6, .5, 10000]
    y_log = np.log10(y)
    p1, success = leastsq(errfunc, p0, args=(x, y_log))
    return p1, success

def logmeans(bins):
    return [np.mean([bins[i], bins[i+1]]) for i in range(len(bins)-1)]


def log_bin_plot(x, y, factor=1.82, col=None, fig=None, ax=None, xlabel=r'$\bar{\tau}$' + ' (days)', ylabel=r'$\langle O | \bar{\tau} \rangle$', title='Overlap as a Function of Mean Inter-Event time', label=None):
    """
    Used for mean waiting time and weight vs Overlap
    """
    bins = binner.Bins(float, min(x)+1, max(x), 'log', factor)
    bin_means, _, _ = binned_statistic(x, y, bins=bins.bin_limits)
    if not fig:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.semilogx(bins.centers, bin_means, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if col:
        cols, _, _ = binned_statistic(x, col, bins=bins.bin_limits)
        ax.scatter(bins.centers, bin_means, c=cols)

    return fig, ax


def lin_bin_plot(x, y, x_bins=20, xlabel=r'$B$', ylabel=r'$\langle O | B \rangle$', title='Overlap as a function of Burstiness', fig=None, ax=None, label=None, arg='-'):
    """
    Used for burstiness VS overlap
    """
    bins = binner.Bins(float, int(np.floor(min(x))), max(x), 'lin', x_bins)
    bin_means, _, _ = binned_statistic(x, y, bins=bins.bin_limits)
    bin_cts, _, _ = binned_statistic(x, y, bins=bins.bin_limits, statistic='count')
    bin_means[bin_cts < 5] = np.nan
    if not fig:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot(bins.centers, bin_means, arg, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig, ax

def plot_overlap_weight_difference(weight_1, overlap_1, weight_2, overlap_2):
    fig, ax = plots.log_bin_plot(weight_1, overlap_2, xlabel=r'$W$', ylabel=r'$\langle O | W \rangle$', title="Overlap as a fucnction of Number \n of calls", label='Company')
    fig, ax = plots.log_bin_plot(weight_2, overlap_2, fig=fig, ax=ax, label='All')


def plot_overlap_difference(simple, extended, bins=35):
    """
    used for plotting overlap calculated with the company set (set 2) vs company plus other users (set 5)
    """
    fig, ax = lin_bin_plot(simple, extended, bins, ylabel=r'$\langle O_{all} | O_{company} \rangle$', title='Overlap for all users VS overlap \n of company users', xlabel=r'$ O_{company}$')
    ax.plot(simple, extended, '.c', alpha=0.03)
    return fig, ax



def loglogheatmap(x, y, z, factor_x=1.5, factor_y=1.45, stat='mean', xlabel=r'$w$ (calls)', ylabel=r'$\bar{\tau}$ (days)', title='Overlap as a function of and Number of Calls and Inter-event Time\n' + r'$\langle O | w, \bar{\tau} \rangle$'):
    x = x+1
    y = y+1
    bins_x = binner.Bins(float, min(x), max(x), 'log', factor_x)
    bins_y = binner.Bins(float, min(y), max(y), 'log', factor_y)
    bin_means, _, _, _ = binned_statistic_2d(x, y, z, statistic=stat, bins=[bins_x.bin_limits, bins_y.bin_limits])
    print(bin_means.shape)
    bin_means = np.nan_to_num(bin_means.T)
    extent = [bins_x.bin_limits[0], bins_x.bin_limits[-1], bins_y.bin_limits[0], bins_y.bin_limits[-1]]
    fig, ax = plt.subplots(1)
    cax = ax.imshow(bin_means, origin="lower") #extent=extent
    #x_ticks = np.linspace(bins_x.bin_limits[0], bins_x.bin_limits[-1],
    #        len(bins_x.bin_limits))
    #ax.set_xticks(x_ticks)
    ax.set_xticklabels(np.int16(bins_x.bin_limits))
    #y_ticks = np.linspace(bins_y.bin_limits[0], bins_y.bin_limits[-1],
    #        len(bins_y.bin_limits))
    #ax.set_yticks(y_ticks)
    ax.set_yticklabels(np.round(bins_y.bin_limits, 1))
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.colorbar(cax)

    return fig, ax


def loglinheatmap(x, y, z, factor_x=1.5, n_bins_y=20, stat='mean', xlabel=r'$w$ (calls)', ylabel=r'$B$', title='Overlap as a function of Number of Calls and Burstiness\n' + r'$\langle O | w, B \rangle$', exp_f=1):
    """
    exp_f contains an "expansion factor". The linear part is ploted [0-1], but we use this expansion factor to depict different things (for instance, =25 for 24 hours)
    """

    y = np.array(y)
    y_alt = max(x)*y/exp_f
    bins_x = binner.Bins(float, min(x), max(x), 'log', factor_x)
    bins_y = binner.Bins(float, min(y_alt), max(y_alt), 'lin', n_bins_y)
    bin_means, _, _, _ = binned_statistic_2d(x, y_alt, z, statistic=stat, bins=[bins_x.bin_limits, bins_y.bin_limits])
    bin_means = np.nan_to_num(bin_means.T)
    extent = [bins_x.bin_limits[0], bins_x.bin_limits[-1], bins_y.bin_limits[0], bins_y.bin_limits[-1]]
    fig, ax = plt.subplots(1)
    cax = ax.imshow(bin_means, extent=extent, origin="lower")
    x_ticks = np.linspace(bins_x.bin_limits[0], bins_x.bin_limits[-1],
            len(bins_x.bin_limits))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(np.int16(bins_x.bin_limits))
    y_ticks = np.linspace(bins_y.bin_limits[0], bins_y.bin_limits[-1],
            len(bins_y.bin_limits))
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([round(yb, 2) for yb in bins_y.bin_limits])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.colorbar(cax)

    return fig, ax


def loglinsbmap(x, y, z, x_factor=1.2, y_bins=40, stat='mean', fig=None, ax=None):
    """
    Heatmap using seaborn

    """
    bins_x = binner.Bins(float, min(x), max(x), 'log', x_factor)
    bins_y = binner.Bins(float, min(y), max(y), 'lin', y_bins)
    bin_means, x_edge, y_edge, _ = binned_statistic_2d(x, y, z, statistic=stat, bins=[bins_x.bin_limits, bins_y.bin_limits])
    bin_means = np.nan_to_num(bin_means.T)
    ker = np.ones((5, 5)); ker[3, 3] = 10.; ker = ker/34.
    bin_means = convolve(bin_means, ker) #Add weights
    bin_means = np.flip(bin_means, 0)
    if fig is not None:
        fig, ax = plt.subplots(1)
    y_ticks = np.flip([str(round(a, 2)) for a in y_edge], 0)
    x_ticks = np.array([str(int(a)) for a in x_edge])
    y_ticks[1:y_bins:2] = ''
    x_ticks[1:len(x_edge):2] = ''
    sb.heatmap(bin_means, ax=ax, xticklabels=x_ticks, yticklabels=y_ticks, robust=True)

    return fig, ax


def paper_heatmaps(x1, x2, x3, y, z, x1_bins=40, x2_bins=30, x3_bins=30, y_bins=40, stat='mean', cbar_kws={}):
    fig, axn = plt.subplots(1, 3, sharey=True)
    cbar_ax = fig.add_axes([.9, .2, .015, .7])

    bins_x1 = binner.Bins(float, min(x1), max(x1), 'lin', x1_bins)
    bins_x2 = binner.Bins(float, min(x2), max(x2), 'lin', x2_bins)
    bins_x3 = binner.Bins(float, min(x3), max(x3), 'lin', x3_bins)
    bins_y = binner.Bins(float, min(y), max(y), 'lin', y_bins)

    bin_means, x_edge, y_edge, _ = binned_statistic_2d(x1, y, z, statistic=stat, bins=[bins_x1.bin_limits, bins_y.bin_limits])
    bin_means, x_edge, y_edge = remove_empty_bins(bin_means.T, x_edge, y_edge)
    bin_means = np.nan_to_num(bin_means)
    ker = np.ones((5, 5)); ker[3, 3] = 10.; ker = ker/34.
    bin_means = convolve(bin_means, ker) #Add weights
    bin_means = np.flip(bin_means, 0)
    y_ticks = np.flip([str(round(a, 2)) for a in y_edge], 0)
    x_ticks_h = np.array([str(round(a, 1)) for a in x_edge])
    xl = len(x_edge)
    x_ticks = ['']*xl
    x_ticks[1:xl:4] = x_ticks_h[1:xl:4]
    y_ticks[1:y_bins:2] = ''
    y_ticks[1:y_bins:3] = ''
    g = sb.heatmap(bin_means, xticklabels=x_ticks, yticklabels=y_ticks, robust=True, square=True, ax=axn[0], cbar=False, vmin=0, vmax=.12)
    g.set(xlabel='(a)    ' + r'$Rank(N^E_{ij})$', ylabel=r'$Rank(w_{ij})$')
    #g.set_tickxlabels(rotation=40)

    bin_means, x_edge, y_edge, _ = binned_statistic_2d(x2, y, z, statistic=stat, bins=[bins_x2.bin_limits, bins_y.bin_limits])
    bin_means, x_edge, y_edge = remove_empty_bins(bin_means.T, x_edge, y_edge)
    bin_means = np.nan_to_num(bin_means)
    ker = np.ones((5, 5)); ker[3, 3] = 10.; ker = ker/34.
    bin_means = convolve(bin_means, ker) #Add weights
    bin_means = np.flip(bin_means, 0)
    y_ticks = np.flip([str(round(a, 2)) for a in y_edge], 0)
    x_ticks_h = np.array([str(round(a, 1)) for a in x_edge])
    xl = len(x_edge)
    x_ticks = ['']*xl
    x_ticks[1:xl:4] = x_ticks_h[1:xl:4]
    y_ticks[1:y_bins:2] = ''
    y_ticks[1:y_bins:3] = ''
    g = sb.heatmap(bin_means, xticklabels=x_ticks, yticklabels=y_ticks, robust=True, square=True, ax=axn[1], cbar=False, vmin=0, vmax=.12)
    g.set(xlabel='(b)    ' + r'$Rank(JSD_{ij})$')

    bin_means, x_edge, y_edge, _ = binned_statistic_2d(x3, y, z, statistic=stat, bins=[bins_x3.bin_limits, bins_y.bin_limits])
    bin_means, x_edge, y_edge = remove_empty_bins(bin_means.T, x_edge, y_edge)
    bin_means = np.nan_to_num(bin_means)
    ker = np.ones((5, 5)); ker[3, 3] = 10.; ker = ker/34.
    bin_means = convolve(bin_means, ker) #Add weights
    bin_means = np.flip(bin_means, 0)
    y_ticks_h = np.flip([str(round(a, 2)) for a in y_edge], 0)
    x_ticks_h = np.array([str(round(a, 1)) for a in x_edge])
    xl, yl = len(x_edge), len(y_edge)
    x_ticks = ['']*xl
    y_ticks = ['']*yl
    x_ticks[1:xl:4] = x_ticks_h[1:xl:4]
    y_ticks[1:yl:3] = y_ticks_h[1:yl:3]
    g = sb.heatmap(bin_means, xticklabels=x_ticks, yticklabels=y_ticks, robust=True, square=True, ax=axn[2], cbar=True, vmin=0, vmax=.12, cbar_ax=cbar_ax, cbar_kws={'label': r'$O_{ij}$'})
    g.set(xlabel='(c)    ' + r'$Rank(B_{ij})$')
    fig.tight_layout(rect=[0, .03, .9, .98])
    fig.savefig('full_run/figs/abstract1.pdf')

    return g


def linlinsbmap(x, y, z, x_bins=40, y_bins=40, stat='mean', z_label=r'$O_{ij}$'):
    """
    Heatmap using seaborn
    """

    bins_x = binner.Bins(float, min(x), max(x), 'lin', x_bins)
    bins_y = binner.Bins(float, min(y), max(y), 'lin', y_bins)
    bin_means, x_edge, y_edge, _ = binned_statistic_2d(x, y, z, statistic=stat, bins=[bins_x.bin_limits, bins_y.bin_limits])

    bin_means, x_edge, y_edge = remove_empty_bins(bin_means.T, x_edge, y_edge)

    bin_means = np.nan_to_num(bin_means)
    ker = np.ones((5, 5)); ker[3, 3] = 10.; ker = ker/34.
    bin_means = convolve(bin_means, ker) #Add weights
    bin_means = np.flip(bin_means, 0)

    y_ticks = np.flip([str(round(a, 2)) for a in y_edge], 0)
    x_ticks = np.array([str(round(a, 2)) for a in x_edge])
    x_ticks[1:x_bins:2] = ''
    y_ticks[1:y_bins:2] = ''
    g = sb.heatmap(bin_means, xticklabels=x_ticks, yticklabels=y_ticks, robust=True, cbar_kws={'label': z_label}, square=True)

    return g


def remove_empty_bins(bin_means, xe, ye):
    xe_n = [xe[0]]
    nan_cols = np.isnan(bin_means).all(0)
    for t, x in zip(nan_cols, xe[1:]):
        if ~t:
            xe_n.append(x)
    bin_means = bin_means[:, ~nan_cols]
    ye_n = [ye[0]]
    nan_rows = np.isnan(bin_means).all(1)
    for t, y in zip(nan_rows, ye[1:]):
        if ~t:
            ye_n.append(y)
    bin_means = bin_means[~nan_rows]
    return bin_means, xe_n, ye_n


def loglogsbmap(x, y, z, x_factor=1.2, y_factor=1.2, stat='mean', fig=None, ax=None):
    """
    Heatmap using seaborn
    """

    bins_x = binner.Bins(float, min(x), max(x), 'log', x_factor)
    bins_y = binner.Bins(float, min(y), max(y), 'log', y_factor)
    bin_means, x_edge, y_edge, _ = binned_statistic_2d(x, y, z, statistic=stat, bins=[bins_x.bin_limits, bins_y.bin_limits])
    bin_means = np.nan_to_num(bin_means.T)
    ker = np.ones((5, 5)); ker[3, 3] = 10.; ker = ker/34.
    bin_means = convolve(bin_means, ker) #Add weights
    bin_means = np.flip(bin_means, 0)
    if fig is not None:
        fig, ax = plt.subplots(1)
    y_ticks = np.flip([str(int(a)) for a in y_edge], 0)
    x_ticks = np.array([str(int(a)) for a in x_edge])
    x_ticks[1:len(x_edge):2] = ''
    y_ticks[1:len(y_edge):3] = ''
    y_ticks[2:len(y_edge):3] = ''
    sb.heatmap(bin_means, ax=ax, xticklabels=x_ticks, yticklabels=y_ticks, robust=True)

    return fig, ax

def linlinheatmap(x, y, z, n_bins_x=30, n_bins_y=30, stat='mean', xlabel=r'$w$ (calls)', ylabel=r'$B$', title='Overlap as a function of Number of Calls and Burstiness\n' + r'$\langle O | w, B \rangle$', exp_f=1):
    """
    exp_f contains an "expansion factor". The linear part is ploted [0-1], but we use this expansion factor to depict different things (for instance, =25 for 24 hours)
    """

    bins_x = binner.Bins(float, min(x), max(x), 'lin', n_bins_x)
    bins_y = binner.Bins(float, min(y), max(y), 'lin', n_bins_y)
    bin_means, _, _ = binned_statistic_2d(x, y, z, statistic=stat, bins=[bins_x.bin_limits, bins_y.bin_limits])
    bin_means = np.nan_to_num(bin_means.T)
    extent = [bins_x.bin_limits[0], bins_x.bin_limits[-1], bins_y.bin_limits[0], bins_y.bin_limits[-1]]
    fig, ax = plt.subplots(1)
    cax = ax.imshow(bin_means, extent=extent, origin="lower")
    x_ticks = np.linspace(bins_x.bin_limits[0], bins_x.bin_limits[-1],
            len(bins_x.bin_limits))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([round(xb, 2) for xb in bins_x.bin_limits])
    y_ticks = np.linspace(bins_y.bin_limits[0], bins_y.bin_limits[-1],
            len(bins_y.bin_limits))
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([round(yb, 2) for yb in bins_y.bin_limits])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.colorbar(cax)
    return fig, ax


def loglinjointdistr(x, y, bins=25, kind='hex', xlim=(1, 2000), ylim=(0, .1), gridsize=(30, 55), xlabel=r'$w_{ij}$', ylabel=r'$O_{ij}$', height=5):
    log_bins = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), bins + 1)
    lin_bins = np.linspace(ylim[0], ylim[1], bins + 1)
    g = sb.jointplot(x, y, kind=kind, xscale='log', marginal_kws=dict(color='w'), xlim=xlim, ylim=ylim, gridsize=gridsize, height=height)

    counts_x = g.ax_marg_x.hist(x, bins=log_bins)
    g.ax_marg_x.set(ylim=(0, max(counts_x[0])))

    counts_y = g.ax_marg_y.hist(y, bins=lin_bins, orientation='horizontal')
    g.ax_marg_y.set(xlim=(0, max(counts_y[0])))
    bin_means, _, _ = binned_statistic(x, y, bins=log_bins)
    log_bins = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), bins)
    g.ax_joint.plot(log_bins, bin_means, 'gray')
    g.ax_joint.set(xlabel=xlabel, ylabel=ylabel)
    return g


def linlinjointdistr(x, y, bins=25, kind='hex', xlim=(-1, 1), ylim=(0, .3), gridsize=(30, 55), xlabel=None, ylabel=r'$O_{ij}$', height=5):
    lin_bins_x = np.linspace(xlim[0], xlim[1], bins + 1)
    lin_bins = np.linspace(ylim[0], ylim[1], bins + 1)
    g = sb.jointplot(x, y, kind=kind, marginal_kws=dict(color='w'), xlim=xlim, ylim=ylim, gridsize=gridsize, height=height)

    counts_x = g.ax_marg_x.hist(x, bins=lin_bins_x)
    g.ax_marg_x.set(ylim=(0, max(counts_x[0])))

    counts_y = g.ax_marg_y.hist(y, bins=lin_bins, orientation='horizontal')
    g.ax_marg_y.set(xlim=(0, max(counts_y[0])))
    bin_means, _, _ = binned_statistic(x, y, bins=lin_bins_x)
    lin_bins_x = np.linspace(xlim[0], xlim[1], bins)
    g.ax_joint.plot(lin_bins_x, bin_means, 'gray')
    g.ax_joint.set(xlabel=xlabel, ylabel=ylabel)
    return g


def loglogjointdistr(x, y, bins=25, kind='hex', xlim=(.01, 1), ylim=(.01, .1), gridsize=(55, 33), xlabel=None, ylabel=r'O_{ij}$', height=5):
    log_bins_x = np.logspace(xlim[0], xlim[1], bins + 1)
    log_bins = np.logspace(ylim[0], ylim[1], bins + 1)
    g = sb.jointplot(x, y, kind=kind, xscale='log', yscale='log', marginal_kws=dict(color='w'), xlim=xlim, ylim=ylim, gridsize=gridsize, height=height)

    counts_x = g.ax_marg_x.hist(x, bins=log_bins_x)
    g.ax_marg_x.set(ylim=(0, max(counts_x[0])))

    counts_y = g.ax_marg_y.hist(y, bins=log_bins, orientation='horizontal')
    g.ax_marg_y.set(xlim=(0, max(counts_y[0])))
    #bin_means, _, _ = binned_statistic(x, y, bins=log_bins_x)
    #log_bins_x = np.logspace(xlim[0], xlim[1], bins)
    #.ax_joint.plot(log_bins_x, bin_means, 'gray')
    g.ax_joint.set(xlabel=xlabel, ylabel=ylabel)
    return g


def cumulative_distribution(x, y, label='', size=100, fig=None, ax=None, xlabel='', ylabel='', title='', arg='.'):
    if not fig:
        fig, ax = plt.subplots(1)
    p_cum, y_mean, sum_x = [], [], 0
    n = float(len(x))
    x = rankdata(x)/n
    x_vals = np.linspace(0, 1 + 1/size, size)
    for v0, v1 in zip(x_vals[:-1], x_vals[1:]):
        ind_x = x < v1
        ind_y = (x >= v0) & ind_x
        y_h = y[ind_y]
        sum_x += len(y_h)
        p_cum.append(sum_x/n)
        y_mean.append(np.mean(y_h))
    ax.plot(p_cum, y_mean, arg, label=label)
    return fig, ax

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 11, # fontsize for x and y labels (was 10)
              'axes.titlesize': 11,
              'text.fontsize': 10, # was 10
              'legend.fontsize': 10, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)


def format_axes(ax):
    SPINE_COLOR = 'gray'
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax


