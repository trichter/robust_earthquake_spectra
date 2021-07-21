# Copyright 2021 Tom Eulenfeld, MIT license

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pickle
from qopen.core import get_pair, Gsmooth
from qopen.rt import G as G_func


def set_gridlabels(ax, i, n, N, xlabel='frequency (Hz)', ylabel=None):
    if i % n != 0 and ylabel:
        plt.setp(ax.get_yticklabels(), visible=False)
    elif i // n == (n - 1) // 2 and ylabel:
        ax.set_ylabel(ylabel)
    if i < N - n and xlabel:
        plt.setp(ax.get_xticklabels(), visible=False)
    elif i % n == (n - 1) // 2 and i >= N - n - 1 and xlabel:
        ax.set_xlabel(xlabel)


def _get_times(tr):
    t0 = tr.stats.starttime - tr.stats.origintime
    return np.arange(len(tr)) * tr.stats.delta + t0


def plot_fits(energies, g0, b, W, R, v0, info, smooth=None,
              smooth_window='bartlett'):
    fs = 250 / 25.4
    plt.figure(figsize=(fs, 0.6*fs))
    tcoda, tbulk, Ecoda, Ebulk, Gcoda, Gbulk = info
    N = len(energies)
    nx, ny = 3, 3
    gs = gridspec.GridSpec(ny, nx, wspace=0.06, hspace=0.08)
    share = None
    if b is None:
        b = 0
    c1 = 'mediumblue'
    c2 = 'darkred'
    c1l = '#8181CD'
    c2l = '#8B6969'
    for i, energy in enumerate(energies):
        evid, station = get_pair(energy)
        ax = plt.subplot(gs[i // nx, i % nx], sharex=share, sharey=share)
        plot = ax.semilogy

        def get_Emod(G, t):
            return R[station] * W[evid] * G * np.exp(-b * t)
        st = energy.stats
        r = st.distance
        t = _get_times(energy) + r / v0 - (st.sonset - st.origintime)

        if smooth:
            plot(t, energy.data_unsmoothed, color='0.7')
        plot(t, energy.data, color=c1l)
        G_ = Gsmooth(G_func, r, t, v0, g0, smooth=smooth,
                     smooth_window=smooth_window)
        Emod = get_Emod(G_, t)
        index = np.argwhere(Emod < 1e-30)[-1]
        Emod[index] = 1e-30

        plot(t, Emod, color=c2l)

        plot(tcoda[i], Ecoda[i], color=c1)
        Emodcoda = get_Emod(Gcoda[i], tcoda[i])
        plot(tcoda[i], Emodcoda, color=c2)

        if tbulk and len(tbulk) > 0:
            plot(tbulk[i], Ebulk[i], 'o', color=c1, mec=c1, ms=4)
            Emodbulk = get_Emod(Gbulk[i], tbulk[i])
            plot(tbulk[i], Emodbulk, 'o', ms=3,
                 color=c2, mec=c2)

        l = '%s\n%dkm' % (station, r / 1000)
        ax.annotate(l, (1, 1), (-5, -5), 'axes fraction',
                    'offset points', ha='right', va='top', size='x-small')

        ylabel = 'spectral energy density $E$ (Jm$^{-3}$Hz$^{-1}$)'
        set_gridlabels(ax, i, nx, N, xlabel='time (s)', ylabel=ylabel)
        kw = dict(color='darkgreen', alpha=0.5, lw=0, zorder=10000)
        ax.axvspan(tcoda[i][0]-4, tcoda[i][0]-0.3, 0.05, 0.08, **kw)
        ax.axvspan(tcoda[i][0]+0.3, tcoda[i][-1], 0.05, 0.08, **kw)
        if share is None:
            share = ax
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_yticks(10. ** np.arange(-11, -5, 2))
    ax.set_xlim((-2, 62))
    ax.set_ylim((1e-13 / 1.5, 1e-6 * 1.5))


if __name__ == '__main__':
    fname = '../qopen/01_go/fits_20186784_04.00Hz-08.00Hz.pkl'
    with open(fname, 'rb') as f:
        tup = pickle.load(f)
    plot_fits(*tup)
    plt.savefig('../figs/qopen_fits_20186784_4-8Hz.pdf', bbox_inches='tight')
