# Copyright 2021 Tom Eulenfeld, MIT license

import matplotlib.pyplot as plt
import numpy as np
from qopen.core import collect_results
from qopen.util import gerr

from util.events import load_qopen


VOGTL = {
    'f': [1.061, 1.5, 2.121, 3.0, 4.243, 6.0, 8.485, 12.0, 16.971, 24.0, 33.941, 48.0],
    'Qsc': [0.005336, 0.004073, 0.002065, 0.00121, 0.000731, 0.000552, 0.000488, 0.000443, 0.000413, 0.000393, 0.000371, 0.000356],
    'Qi': [0.011312, 0.010018, 0.006624, 0.004049, 0.002507, 0.001673, 0.001213, 0.000956, 0.000761, 0.000605, 0.000484, 0.000388]}

BACHURA2016 = {
    'f': [1.5, 3, 6, 9, 12, 18],
    'Qsc': [0.00341, 0.001949, 0.000957, 0.000622, 0.000485, 0.000294],
    'Qi': [0.005743, 0.00185, 0.001167, 0.000859, 0.000725, 0.000546]}


def Q(q, obs='sc', v0=None, error=False):
    if v0 is None:
        v0 = np.mean([ev['v0'] for ev in q['events'].values()])
        #v0 = q['config']['v0']
        print(v0)
    freq = np.array(q['freq'])
    if obs == 'sc':
        mean = np.array(q['g0']) * v0 / (2 * np.pi * freq)
    else:
        mean = np.array(q['b']) / (2 * np.pi * freq)
    if not error:
        return freq, mean
    q = collect_results(q)
    if obs == 'sc':
        vals = np.array(q['g0']) * v0 / (2 * np.pi * freq)
    else:
        vals = np.array(q['b']) / (2 * np.pi * freq)
    mean2, err1, err2 = gerr(vals, axis=0, robust=True)
    np.testing.assert_allclose(mean, mean2)
    return freq, mean, (err1, err2)


def _vogtland_results2Q():
    q = load_qopen('../data/vogtland_results.json')
    freq, Qsc = Q(q, 'sc')
    _, Qi = Q(q, 'i')
    freq = np.round(freq, 3).tolist()
    Qsc = np.round(Qsc, 6).tolist()
    Qi = np.round(Qi, 6).tolist()
    print(f"""VOGTL = {{
              'f': {freq},
              'Qsc': {Qsc},
              'Qi': {Qi}}}""")


def plotQ():
    q1 = load_qopen('Q')
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharex=ax1, sharey=ax1)


    freq, Qsc = Q(q1, 'sc')
    _, Qi = Q(q1, 'i')
    print(freq)
    print(Qsc)
    print(Qi)
    print(1/(Qi + Qsc))

    ax1.errorbar(*Q(q1, 'sc', error=True), marker='.', ms=6, label='this study', zorder=5)
    ax1.loglog(VOGTL['f'], VOGTL['Qsc'], label='Gaebler et al. 2015')
    ax1.plot(BACHURA2016['f'], BACHURA2016['Qsc'], label='Bachura and Fischer, 2016')
    ax1.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'scattering strength $Q_{\mathrm{sc}}^{-1}$')

    ax2.errorbar(*Q(q1, 'i', error=True), marker='.', ms=6, zorder=5)
    ax2.plot(VOGTL['f'], VOGTL['Qi'])
    ax2.plot(BACHURA2016['f'], BACHURA2016['Qi'])
    ax2.set_xlabel('frequency (Hz)')
    ax2.set_ylabel(r'intrinsic attenuation $Q_{\mathrm{intr}}^{-1}$')
    ax1.legend(loc='upper right', fontsize='small')
    ax1.set_xticks((1, 10, 100))
    ax1.set_xticklabels(('1', '10', '100'))
    ax2.set_xticks((1, 10, 100))
    ax2.set_xticklabels(('1', '10', '100'))
    fig.tight_layout()
    fig.savefig('../figs/Q.pdf')


if __name__ == '__main__':
    plotQ()
