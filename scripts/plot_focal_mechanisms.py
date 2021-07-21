# Copyright 2021 Tom Eulenfeld, MIT license

import matplotlib.pyplot as plt
import numpy as np
from obspy.imaging.beachball import beach

from util.events import get_bounds2, load_grond


def plot_grond_focal_mechanisms(ax=None, color_coded=False):
    gevents = load_grond()
    if ax is None:
        fig = plt.figure(figsize=(2, 6) if color_coded else (2, 2))
        ax = fig.add_subplot(111, aspect=1)
    else:
        fig = None
    ax.axis('off')
    if color_coded:
        bounds, colors, labels = get_bounds2()
    for ev in gevents.values():
        mt = ev['grond'].moment_tensor_dc.m6_up_south_east()
        if color_coded:
            i = np.searchsorted(bounds, ev['mpl'])
            xy = (0, - i * 200)
            color = colors[i]
        else:
            xy = (0, 0)
            color = 'k'
        b = beach(mt, linewidth=1, edgecolor=color, alpha=0.1, xy=xy,
                  width=180, nofill=True, zorder=-100)
        ax.add_collection(b)
    if color_coded:
        for i in range(3):
            ax.annotate(labels[i], (0, 100 - i * 200), size='x-small',
                        ha='right', va='top')
    ax.set_rasterization_zorder(-50)
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100 - color_coded * 400, 100)
    if fig is not None:
        fig.savefig('../figs/focal_mechanisms.pdf', bbox_inches='tight',
                    pad_inches=0.1, dpi=600)


if __name__ == '__main__':
    plot_grond_focal_mechanisms(color_coded=True)
