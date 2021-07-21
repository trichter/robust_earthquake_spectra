# Copyright 2021 Sebastian Heimann, Tom Eulenfeld, MIT license

from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from pyrocko.moment_tensor import kagan_angle, MomentTensor

from util.events import load_grond


def plot_std_vs_mag(events):
    fig = plt.figure()
    ax = None
    mags = []
    stds = defaultdict(list)
    means = defaultdict(list)
    for ev in events.values():
        if 'grond_stats' in ev:
            mags.append(ev['ref'].magnitude)
            for s in ev['grond_stats']:
                if not s.name.startswith('rm'):
                    stds[s.name].append(s.std)
                    means[s.name].append(s.mean)
    for i, prop in enumerate(stds.keys()):
        print(f'median mean for {prop}: {np.median(means[prop]):.3f}')
        print(f'median std for {prop}: {np.median(stds[prop]):.3f}')
        ax = fig.add_subplot(len(stds), 1, 1+i, sharex=ax)
        ax.scatter(mags, stds[prop])
        ax.set_xlabel('magnitude')
        ax.set_ylabel(f'{prop} std')


def calc_kagan(events):
    medkagans = []
    for ev in events.values():
        if 'grond_ensemble' in ev:
            kagans = [kagan_angle(sol.moment_tensor, ev['grond'].moment_tensor)
                      for sol in ev['grond_ensemble']]
            medkagans.append(np.median(kagans))
    print(f'median of median of kagan angles: {np.median(medkagans):.1f}')


if __name__ == '__main__':
    events = load_grond(load_stats=True)
    plot_std_vs_mag(events)
    print('load ensembles this might take a while...')
    events = load_grond(load_ensembles=True)
    calc_kagan(events)
    print(MomentTensor.from_values([268, 59, -167]).both_strike_dip_rake())
