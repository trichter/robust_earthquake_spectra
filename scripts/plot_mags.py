# Copyright 2021 Tom Eulenfeld, MIT license

import matplotlib.pyplot as plt
import numpy as np
from pyrocko.moment_tensor import eigh_check
from qopen.imaging import _secondary_yaxis_seismic_moment
from qopen.source import moment_magnitude
from statsmodels.regression.linear_model import OLS, WLS
from statsmodels.robust.robust_linear_model import RLM

from util.events import get_bounds2, load_qopen_grond_sds


MwQs = r'$M_{\rm{wQ}}$'
MwOs = r'$M_{\rm{wO}}$'
MwGs = r'$M_{\rm{wG}}$'
Mls = r'$M_{\rm{L}}$'
Mws = r'$M_{\rm{w}}$'

OUT = '../figs/'


def linear_fit(y, x, m=None, method='robust', **kw):
    """Linear fit between x and y

    :param y,x: data
    :param m: fix slope at specific value
    :param method: one of ('least_squares', 'weighted', 'robust')
    :return: slope a and intercept b of y = ax + b
    """
    Model = RLM if method == 'robust' else WLS if method == 'weighted' else OLS
    if m is None:
        X = np.empty((len(y), 2))
        X[:, 0] = x
        X[:, 1] = 1
        res = Model(y, X, **kw).fit()
        return res.params, res
    else:
        X = np.ones(len(y))
        res = Model(np.array(y) - m * np.array(x), X, **kw).fit()
        return m, res.params[0]


def residuals(stackh=False):
    res = [(evid,
            ev['qopen']['Mcat'] - moment_magnitude(ev['grond'].moment_tensor.moment),
            ev['qopen']['Mw'] - moment_magnitude(ev['grond'].moment_tensor.moment),
            ev['grond'].north_shift,
            ev['grond'].east_shift,
            ev['grond'].depth,
            ev['grond'].depth - ev['ref'].depth,
            ev['grond'].time - ev['ref'].time,
            ev['grond'].moment_tensor.strike2,
            ev['grond'].moment_tensor.dip2,
            ev['grond'].moment_tensor.rake2,
            ev['grond'].moment_tensor.standard_decomposition()[2],
            ev['mpl']
            )
           for evid, ev in load_qopen_grond_sds().items()]
    evids, M1, M2, north, east, depth, depthshift, time, strike, dip, rake, clvd_decomp, mpl = zip(*res)

    if stackh:
        bounds, colors, _ = get_bounds2()
        mpl = np.array(mpl)
        ind1 = mpl<bounds[0]
        ind2 = np.logical_and(mpl>=bounds[0], mpl<=bounds[1])
        ind3 = mpl>bounds[1]
        def stack_hist(ax, data, bins=None):
            data = np.array(data)
            data = [data[ind1], data[ind2], data[ind3]]
            ax.hist(data, bins=bins, color=colors, histtype='barstacked')
    else:
        def stack_hist(ax, data, bins=None):
            ax.hist(data, bins=bins)

    fig, ax = plt.subplots(2, 4, figsize=(12, 6))
    stack_hist(ax[0, 0], north, bins=np.linspace(-1500, 1500, 26))
    ax[0, 0].set_xlabel('north shift (m)')
    stack_hist(ax[0, 1], east, bins=np.linspace(-1500, 1500, 26))
    ax[0, 1].set_xlabel('east shift (m)')
    #ax[0, 2].hist(depthshift, bins=np.linspace(0, 2500, 26))
    stack_hist(ax[0, 2], depthshift, bins=np.linspace(-1500, 1500, 26))
    ax[0, 2].set_xlabel('depth shift (m)')
    #ax[0, 3].hist(time, bins=np.linspace(-0.5, 0.2, 21))
    stack_hist(ax[0, 3], time, bins=np.linspace(-0.3, 0.3, 21))
    ax[0, 3].set_xlabel('time shift (s)')
    stack_hist(ax[1, 0], strike, bins=np.linspace(0, 360, 37))
    print(f'median strike: {np.median(strike):.0f}')
    ax[1, 0].set_xlabel('strike (°)')
    stack_hist(ax[1, 1], dip, bins=np.linspace(0, 90, 37))
    print(f'median dip: {np.median(dip):.0f}')
    ax[1, 1].set_xlabel('dip (°)')
    stack_hist(ax[1, 2], rake, bins=np.linspace(-180, 180, 37))
    rake2 = np.array(rake)
    rake2[rake2>0] = rake2[rake2>0] - 360
    print(f'median rake: {np.median(rake2):.0f}')
    ax[1, 2].set_xlabel('rake (°)')
    ratios = []
    for res in clvd_decomp:
        # inspired by grond/src/problems/cmt/problem.py
        ratio_clvd, m_clvd = res[1:3]
        evals, evecs = eigh_check(m_clvd)
        ii = np.argmax(np.abs(evals))
        ratio_clvd = ratio_clvd * np.sign(evals[ii])
        ratios.append(ratio_clvd)
    stack_hist(ax[1, 3], ratios, bins=np.linspace(-0.2, 0.2, 26))
    print(f'median clvd ratio: {np.median(ratios):.3f}')
    ax[1, 3].set_xlabel('relative contribution of CLVD\ncomponent to seismic moment')
    for i in range(8):
        axl = ax[i // 4, i % 4]
        label = 'abcdefgh'[i] + ')'
        axl.annotate(label, (0, 1), (5, -5), 'axes fraction', 'offset points',
                     va='top', size='large')
    fig.tight_layout()
    fname = OUT + 'residuals.pdf'
    if stackh:
        fname = fname.replace('.pdf', '_stacked.pdf')
    fig.savefig(fname)


def plot_resid(fig, ax, res, r=0.2):
    bbox = ax.get_position()
    nbbox = [bbox.x0 + 0.7 * bbox.width, bbox.y0 + 0.1 * bbox.height,
             0.2 * bbox.width, 0.2 * bbox.height]
    axn = fig.add_axes(nbbox)
    axn.hist(res.resid, 15, range=(-r, r))
    label = f'RMSE={res.mse_resid**0.5:.3f}'
    ax.annotate(label, (0.6, 0.32), xycoords='axes fraction')
    axn.set_xticks([-r, 0, r])
    axn.set_xticklabels([str(-r), '0', str(r)])
    axn.spines['left'].set_visible(False)
    axn.spines['right'].set_visible(False)
    axn.spines['top'].set_visible(False)
    axn.tick_params(left=False, labelleft=False)


def compare_mags():
    events = [(evid, ev['qopen']['Mcat'], ev['qopen']['Mw'],
               moment_magnitude(ev['grond'].moment_tensor.moment),
               moment_magnitude(ev['M0_onset'])
               )
              for evid, ev in load_qopen_grond_sds().items()]

    evids, Ml, MwQ, MwG, MwO = map(np.array, zip(*events))
    print(evids[Ml>3.5])
    print(Ml[evids=='20186784'])
    method = 'least squares'
    # method = 'robust'
    mag = np.array([1.8, 4])
    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(133, aspect='equal')
    ax2 = fig.add_subplot(131, sharex=ax1, sharey=ax1, aspect='equal')
    ax3 = fig.add_subplot(132, sharex=ax1, sharey=ax1, aspect='equal')
    ax1.plot(MwQ, Ml, 'x', ms=4)#, color='0.3')
    (m, b), res1 = linear_fit(Ml[MwQ>=2.3], MwQ[MwQ>=2.3], method=method)
    label = f'{Mls} = {m:.2f}{MwQs} {b:+.2f}\n{MwQs} = {1/m:.2f}{Mls} {-b/m:+.2f}'
    (m2, b2), res1b = linear_fit(MwQ, Ml, method=method)
    label2 = f'{MwQs} = {m2:.2f}{Mls} {b2:+.2f}'
    print(label)
    print(label2)

    ax1.plot(mag, m * mag + b, 'k', label=label, alpha=0.7)
    ax1.plot(mag, mag / m2 - b2 / m2, '--k', label=label2, alpha=0.5)
    ax1.legend()

    ax2.plot(MwQ, MwG, 'x', ms=4)#, color='0.3')
    (m, b), res2 = linear_fit(MwG, MwQ, method=method )
    _, b2 = linear_fit(MwG, MwQ, m=1, method=method)
    label = f'{MwGs} = {m:.2f}{MwQs} {b:+.2f}'
    label2 = f'{MwGs} = {MwQs} {b2:+.2f}'
    print(label)
    print(label2)
    ax2.plot(mag, mag * m + b, 'k', label=label, alpha=0.7)
    ax2.plot(mag, mag + b2, '--k', label=label2, alpha=0.7)
    ax2.legend()

    ax3.plot(MwQ, MwO, 'x', ms=4)#, color='0.3')
    (m, b), res3 = linear_fit(MwO, MwQ, method=method)
    _, b2 = linear_fit(MwO, MwQ, m=1, method=method)
    label = f'{MwOs} = {m:.2f}{MwQs} {b:+.2f}'
    label2 = f'{MwOs} = {MwQs} {b2:+.2f}'
    print(label)
    print(label2)
    ax3.plot(mag, mag * m + b, 'k', label=label, alpha=0.7)
    ax3.plot(mag, mag + b2, '--k', label=label2, alpha=0.7)
    ax3.legend()
    akwargs = dict(xy=(-0.1, 1.02), xycoords='axes fraction', size='large')
    ax2.annotate('a)', **akwargs)
    ax3.annotate('b)', **akwargs)
    ax1.annotate('c)', **akwargs)
    ax1.set_xlabel(f'Qopen moment magnitude {MwQs}')
    ax2.set_xlabel(f'Qopen moment magnitude {MwQs}')
    ax3.set_xlabel(f'Qopen moment magnitude {MwQs}')
    ax1.set_ylabel(f'WEBNET local magnitude {Mls}')
    ax2.set_ylabel(f'Grond moment magnitude {MwGs}')
    ax3.set_ylabel(f'onset moment mangiutde {MwOs}')
    ax1.set_ylim(ax1.get_xlim())
    fig.tight_layout()
    plot_resid(fig, ax1, res1, r=0.5)
    plot_resid(fig, ax2, res2)
    plot_resid(fig, ax3, res3)
    fig.savefig(OUT + 'mags_comparison.pdf', bbox_inches='tight')


def comparison_Mw_vs_Ml():
    Ml = np.array((1, 4))
    relations = {
            'this study, fig 9c, continous line': (0.95, 0.16),
            'this study, fig 9c, stroken line': (0.84, 0.39),
            'Jakoubková et al., 2018': (0.73, 0.66),
            'Michalek & Fischer, 2013': (0.92, 0.80),
            'Horalek & Silny, 2013': (0.75, 0.45),
            }
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for study, r in relations.items():
        func = lambda Ml: r[0] * Ml + r[1]
        ls = '--' if 'stroken' in study else '-'
        ax.plot(Ml, func(Ml), label=f'{Mws}={r[0]:.2f}{Mls}+{r[1]:.2f}, {study}', ls=ls)
    ax.legend()
    ax.set_xlabel(f'WEBNET local magnitude {Mls}')
    ax.set_ylabel(f'moment magnitude {Mws}')
    _secondary_yaxis_seismic_moment(ax)
    fig.savefig(OUT + 'MwMl.pdf', bbox_inches='tight')


if __name__ == '__main__':
    compare_mags()
    residuals()
    residuals(stackh=True)
    comparison_Mw_vs_Ml()
