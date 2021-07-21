# Copyright 2021 Tom Eulenfeld, MIT license

import json
import matplotlib.pyplot as plt
import numpy as np
from qopen.imaging import (plot_all_sds, plot_sites,
                           _secondary_yaxis_seismic_moment)
from qopen.source import moment_magnitude
from scipy.optimize import minimize
from scipy.stats import linregress

from plot_maps import get_norm
from prepare_data import load_full_catalog


Mws = r'$M_{\rm w}$'
fcs = r'$f_{\rm c}$'

QOPENEVENTRESULTS = '../qopen/03_mag/results.json'
QOPENEVENTRESULTS2 = '../qopen/04_mag_nconst/results.json'
QOPENSITERESULTS = '../qopen/02_sites/results.json'


def _linear_fit_L1(y, x, m0, b0):
    def cost_function(params, x, y):
        m, b = params
        return np.sum(np.abs(y - m * x - b))
    out = minimize(cost_function, (m0, b0), args=(x, y))
    return out.x


def get_cmap():
    allevents = load_full_catalog()
    cmap = plt.get_cmap('turbo')
    norm, _ = get_norm()
    colors = {str(ev.resource_id).split('/')[-1]:
              cmap(norm(ev.origins[0].time.matplotlib_date))
              for ev in allevents}
    return colors


def plot_all_sds_wrapper(plot_fm=False):
    with open(QOPENEVENTRESULTS) as f:
        results = json.load(f)
    colors = get_cmap()
    plot_all_sds(results, nx=8, annotate=False, colors=colors, figsize=(12, 30))
    fig = plt.gcf()
    fig.axes[0].set_yticks((10**10, 10**12, 10**14))
    fig.axes[0].set_xticks((1, 10, 100))
    fig.axes[0].set_xticklabels(('1', '10', '100'))

    if plot_fm:
        from pyrocko.plot import beachball
        from util.events import get_norm, load_grond
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cmap = plt.get_cmap('turbo')
        norm, _ = get_norm()
        gevents = load_grond()
        gevents = sorted(gevents.values(), key= lambda ev: ev['ref'].time)
        for ax, ev in zip(fig.axes, gevents):
            axins = inset_axes(ax, width=0.33, height=0.33, borderpad=0.2)
            axins.axis('off')
            beachball.plot_beachball_mpl(
                ev['grond'].moment_tensor.m6(), axins,
                beachball_type='dc',
                size=10, size_units='data',
                color_t=cmap(norm(ev['mpl'])),
                linewidth=0.5,
                zorder=10
                )
            axins.set_xlim(-5.1, 5.1)
            axins.set_ylim(-5.1, 5.1)
            axins.set_rasterization_zorder(20)
    fig.savefig('../figs/all_sds.pdf', dpi=200)


def fc_vs_Mw_vs_n():
    with open(QOPENEVENTRESULTS) as f:
        results = json.load(f)
    vals = [(evres['Mw'], evres['fc'], evres['n']) for evres in
            results['events'].values() if 'Mw' in evres]
    Mw, fc, n = map(np.array, zip(*vals))
    print(f'high frequency fall-of mean: {np.mean(n):.3f}  '
          f'median: {np.median(n):.3f}  std: {np.std(n):.3f}')
    fig = plt.figure(figsize=(12,5))
    ax1 = fig.add_subplot(131)
    ax1.plot(Mw, fc, 'x')
    ax1.set_xlabel(f'moment magnitude {Mws}')
    ax1.set_ylabel(f'corner frequency {fcs} (Hz)')
    ax1.set_yscale('log')
    ax1.set_yticks([3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ax1.set_yticklabels(['3', '', '5', '', '7', '', '', '10', '', '12'])
    ax2 = fig.add_subplot(132)
    ax2.hist(n, bins=np.arange(1.525, 2.26, 0.05), rwidth=0.9)
    ax2.set_xlabel('high frequency fall-off $n$')
    ax3 = fig.add_subplot(133)
    ax3.plot(Mw, n, 'x')
    ax3.set_xlabel(f'moment magnitude {Mws}')
    ax3.set_ylabel('high frequency fall-off $n$')
    plt.tight_layout()
    fig.savefig('../figs/fc_vs_Mw_vs_n.pdf', bbox_inches='tight', pad_inches=0.1)


def fc2Mw(stress_drop, fc):
    r = 3500 * 0.21 / np.array(fc)  # Madariaga (1976) for S waves
    M0 = 16 * r ** 3 * stress_drop / 7
    return moment_magnitude(M0)


def fc2stress_drop(fc, M0):
    r = 3500 * 0.21 / np.array(fc)
    stress_drop = 7 * M0 / (16 * r ** 3)
    return stress_drop


def fc_vs_Mw_fixed_n(show_n=True, color_coded=False, color_linreg=False):
    with open(QOPENEVENTRESULTS2) as f:
        results = json.load(f)
    vals = [(evres['Mw'], evres['fc']) for evres in
            results['events'].values() if 'Mw' in evres]
    Mw, fc = map(np.array, zip(*vals))
    fig = plt.figure(figsize=(10 if show_n else 5, 5))
    ax1 = fig.add_axes([0.45, 0.1, 0.45, 0.85]) if show_n else fig.add_subplot(111)
    if color_coded or color_linreg:
        from util.events import get_norm, get_bounds2, load_grond
        cmap = plt.get_cmap('turbo')
        norm, _ = get_norm()
        bounds, colors, labels = get_bounds2()
        gevents = load_grond()
        mpl = [gevents[evid]['mpl'] for evid, evres in
                results['events'].items() if 'Mw' in evres]
        mpl = np.array(mpl)
        ind1 = mpl < bounds[0]
        ind2 = np.logical_and(mpl >= bounds[0], mpl <= bounds[1])
        ind3 = mpl > bounds[1]
    if not color_coded:
        ax1.plot(fc, Mw, 'x', color='k' if color_linreg else 'C0')
    else:
        if color_coded == 'mean':
            ax1.scatter(fc[ind1], Mw[ind1], marker='x', color=colors[0])
            ax1.scatter(fc[ind2], Mw[ind2], marker='x', color=colors[1])
            ax1.scatter(fc[ind3], Mw[ind3], marker='x', color=colors[2])
        else:
            ax1.scatter(fc, Mw, marker='x', c=mpl, cmap=cmap, norm=norm)
    m, b, _, _, m_stderr = linregress(np.log10(fc), Mw)
    print(f'L2 fit, fc independent variable: M0 ∝ fc^{1.5*m:.2f}+-{1.5*m_stderr:.2f}')
    m, b = _linear_fit_L1(Mw, np.log10(fc), m, b)
    print(f'L1 fit, fc independent variable: M0 ∝ fc^{1.5*m:.2f}')
    m2, b2 = _linear_fit_L1(np.log10(fc), Mw, m, b)
    m, b = 1 / m2, -b2 / m2
    print(f'L1 fit, Mw independent variable: M0 ∝ fc^{1.5*m:.2f}')
    m3, b3, _, _, m3_stderr = linregress(Mw, np.log10(fc))
    m, b = 1 / m3, -b3 / m3
    m_stderr = m3_stderr / m3 ** 2
    print(f'L2 fit, Mw independent variable: M0 ∝ fc^{1.5*m:.2f}+-{1.5*m_stderr:.2f}')
    fclim = np.array([2, 5, 8, 15])
    label = f'$M_0$ ∝ {fcs}$^{{{1.5*m:.2f} \pm {1.5*m_stderr:.2f}}}$'
    ax1.plot(fclim, np.log10(fclim)*m+b, zorder=-1, label=label,
             color='k' if color_coded or color_linreg else 'C1')
    if color_linreg:
        for i in range(3):
            ind = [ind1, ind2, ind3][i]
            m2, b2, _, _, m2_stderr = linregress(Mw[ind], np.log10(fc[ind]))
            m, b = 1 / m2, -b2 / m2
            m_stderr = m2_stderr / m2 ** 2
            print(f'L2 fit, Mw independent variable, {labels[i]}: M0 ∝ fc^{1.5*m:.2f}+-{1.5*m_stderr:.2f}')
            label = f'$M_0$ ∝ {fcs}$^{{{1.5*m:.2f} \pm {1.5*m_stderr:.2f}}}$'
            ax1.plot(fclim, np.log10(fclim)*m+b, zorder=-1, label=label,
                     color=colors[i], alpha=0.5 if color_coded else 0.8)
    kw = dict(ls='--', color='0.5', zorder=-1)
    ax1.plot(fclim, fc2Mw(0.1e6, fclim), **kw)
    ax1.plot(fclim, fc2Mw(1e6, fclim), label=f'$M_0$ ∝ {fcs}$^{{-3}}$', **kw)
    ax1.plot(fclim, fc2Mw(10e6, fclim), **kw)
    kw = dict(rotation=-35, rotation_mode='anchor',
              xytext=(0, -15), textcoords='offset points')
    ax1.annotate('0.1 MPa', (2.4, fc2Mw(0.1e6, 2.4)), **kw)
    ax1.annotate('1 MPa', (2.4, fc2Mw(1e6, 2.4)), **kw)
    ax1.annotate('10 MPa', (11, fc2Mw(10e6, 11)), **kw)
    fc1 = 10 ** (moment_magnitude(1e13) / m - b / m)
    sd1 = fc2stress_drop(fc1, 1e13) / 1e6
    print(f'fc={fc1:.1f}Hz and stress drop {sd1:.1f} Mpa for M0=1e13 Nm')
    ax1.set_ylim(1.7, None)
    ax1.set_xlim(2.25, 14)
    ax1.set_ylabel(f'moment magnitude {Mws}')
    ax1.set_xlabel(f'corner frequency {fcs} (Hz)')
    ax1.set_xscale('log')
    ax1.set_xticks([3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ax1.set_xticklabels(['3', '', '5', '', '7', '', '', '10', '', '12'])
    ax1.legend()
    _secondary_yaxis_seismic_moment(ax1)
    if show_n:
        with open(QOPENEVENTRESULTS) as f:
            results = json.load(f)
        vals = [(evres['Mw'], evres['fc'], evres['n']) for evres in
                results['events'].values() if 'Mw' in evres]
        Mw, fc, n = map(np.array, zip(*vals))
        print(f'high frequency fall-of mean: {np.mean(n):.3f}  '
              f'median: {np.median(n):.3f}  std: {np.std(n):.3f}')
        ax3 = fig.add_axes([0.05, 0.1, 0.3, 0.85])
        ax3.hist(n, bins=np.arange(1.525, 2.26, 0.05), rwidth=0.9)
        ax3.set_xlabel('high frequency fall-off $n$')
        akwargs = dict(xy=(0.02, 0.95), xycoords='axes fraction', size='large')
        ax1.annotate('b)', **akwargs)
        ax3.annotate('a)', **akwargs)
        ax1.set_ylim(None, 3.9)
    fname = '../figs/n_and_fc_vs_Mw.pdf' if show_n else '../figs/fc_vs_Mw.pdf'
    if color_coded:
        fname = fname.replace('.pdf', '_color.pdf')
    if color_coded == 'mean':
        fname = fname.replace('.pdf', '2.pdf')
    if color_linreg:
        fname = fname.replace('.pdf', '_clinreg.pdf')
    fig.savefig(fname, bbox_inches='tight', pad_inches=0.1)


def plot_sites_():
    with open(QOPENSITERESULTS) as f:
        results = json.load(f)
    plot_sites(results, nx=5, figsize=(12, 5), ylim=(10**-1.2, 10**1.2),
               ylabel=None)
    fig = plt.gcf()
    fig.axes[0].set_xticks((1, 10, 100))
    fig.axes[0].set_xticklabels(('1', '10', '100'))
    for ax in fig.axes[:-2]:
        ax.axhspan(0.5, 2, color='0.85', zorder=-9)
        ax.axhline(1, color='0.75', zorder=-5)
    fig.supylabel('energy site amplification', x=0.07)
    fig.savefig('../figs/sites.pdf')


if __name__ == '__main__':
    plot_all_sds_wrapper(plot_fm=True)
    fc_vs_Mw_vs_n()
    fc_vs_Mw_fixed_n(show_n=False)
    fc_vs_Mw_fixed_n(color_coded='mean', color_linreg=True)
    fc_vs_Mw_fixed_n(color_linreg=True)
    plot_sites_()
    # fc_vs_Mw_fixed_n()
    # fc_vs_Mw_fixed_n(color_coded=True)
    # fc_vs_Mw_fixed_n(color_coded=True, color_linreg=True)
    sm1 = moment_magnitude(0.6, inverse=True)
    sm2 = moment_magnitude(0, inverse=True)
    fac = round(sm2 ** (-1/4.7) / sm1 ** (1/-4.7), 2)
    print(f'offset of 0.6 in magnitude -> factor {fac:.2f} in corner frequency')
