# Copyright 2021 Tom Eulenfeld, MIT license

from matplotlib.dates import DateFormatter
from matplotlib.colorbar import ColorbarBase
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy import UTCDateTime as UTC
from pyrocko.plot import beachball
from qopen.source import moment_magnitude

from prepare_data import load_full_catalog
from util.events import events2lists, load_grond, load_qopen, get_norm
from util.imaging import convert_coords2km


LATLON0 = (50.25, 12.45)


def Ml2cumM0(evids, mags):
    res = load_qopen('events')
    M0s = []
    for evid, mag in zip(evids, mags):
        try:
            M0 = res['events'][evid]['M0']
        except KeyError:
            Mw = 0.95 * mag + 0.16
            M0 = moment_magnitude(Mw, inverse=True)
        M0s.append(M0)
    return np.cumsum(M0s)


def plot_events_stations_map_depth(events, inv=None, figsize=(8,8), out=None, show=True, cmap='turbo',
                                   dpi=300, convert_coords=False, all_events=None, colorbar=True, label=True, ms=9,
                                   plot_fm=False, plot_color_ax=False):
    """
    Modified from github.com/trichter/inter_source_interferometry
    """
    id_, time, lon, lat, dep, mag, *_ = events
    mag = np.array(mag)
    mpl = np.array([t.matplotlib_date for t in time])
    if all_events:
        id2, time2, lon2, lat2, dep2, mag2, *_ = all_events
        mag2 = np.array(mag2)
        mpl2 = np.array([t.matplotlib_date for t in time2])
        lon2 = np.array(lon2)[mag2<1.9]
        lat2 = np.array(lat2)[mag2<1.9]
        dep2 = np.array(dep2)[mag2<1.9]
        mag3 = mag2[mag2<1.9]
        mpl3 = mpl2[mag2<1.9]
    fig = plt.figure(figsize=figsize)
    if plot_fm:
        ax3 = fig.add_axes((0.1, 0.2, 0.4, 0.7))
    else:
        ax3 = fig.add_axes((0.1, 0.5, 0.3, 0.4))
    if not plot_fm:
        ax4 = fig.add_axes((0.42, 0.5, 0.35, 0.4), sharey=ax3)
        ax5 = fig.add_axes((0.1, 0.5-0.37, 0.3, 0.35), sharex=ax3)
        def _on_lims_changed(ax, boo=[True]):
            if boo[0]:
                boo[0] = False
                if ax == ax5:
                    ax4.set_xlim(ax5.get_ylim()[::-1])
                if ax == ax4:
                    ax5.set_ylim(ax4.get_xlim()[::-1])
                boo[0] = True

    if not plot_fm:
        ax5.invert_yaxis()
        ax5.callbacks.connect('ylim_changed', _on_lims_changed)
        ax4.callbacks.connect('xlim_changed', _on_lims_changed)
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position("right")
        ax3.xaxis.tick_top()
        ax3.xaxis.set_label_position("top")
        ax4.xaxis.tick_top()
        ax4.xaxis.set_label_position("top")


    if convert_coords:
        latlon0 = None if convert_coords==True else convert_coords
        x, y = zip(*convert_coords2km(list(zip(lat, lon)), latlon0=latlon0))
        if all_events:
            x2, y2 = zip(*convert_coords2km(list(zip(lat2, lon2)), latlon0=latlon0))
    else:
        x, y = lon, lat
        if all_events:
            x2, y2 = lon2, lat2

    norm, bounds = get_norm(extend=True)
    cmap = plt.get_cmap(cmap)
    if all_events:
        ax3.scatter(x2, y2, 4, color='0.6')
        ax4.scatter(dep2, y2, 4, color='0.6')
        ax5.scatter(x2, dep2, 4, color='0.6')
    if plot_fm:
        for xx, yy, dd, t, mt in zip(x, y, dep, mpl, plot_fm):
            beachball.plot_beachball_mpl(
                mt.m6(), ax3,
                beachball_type='dc',
                size=10,
                edgecolor='black',
                position=(xx, yy),
                color_t=cmap(norm(t)),
                linewidth=0.5,
                zorder=-100)
    else:
        ax3.scatter(x, y, ms, mpl, cmap=cmap, norm=norm)
        ax4.scatter(dep, y, ms, mpl, cmap=cmap, norm=norm)
        ax5.scatter(x, dep, ms, mpl, cmap=cmap, norm=norm)
    if colorbar:
        if plot_fm:
            ax7 = fig.add_axes([0.1, 0.08, 0.34, 0.02])
        else:
            ax7 = fig.add_axes([0.56, 0.42, 0.34, 0.02])
        cbar = ColorbarBase(ax7, cmap=cmap, norm=norm, orientation='horizontal', format=DateFormatter('%Y-%m-%d'), extend='max')
        cbar.ax.xaxis.set_minor_locator(MultipleLocator(1))
        cbar.set_ticks([UTC(f'2018-05-{day}').matplotlib_date for day in (10, 13, 16, 19, 22, 25, 26)])
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels()[:-1] + ['until 2018-06-19'], rotation=60, ha='right', rotation_mode='anchor')
    if plot_color_ax:
        figy0 = 0.18
        axc = fig.add_axes((0.5, figy0, 0.3, 0.25))
        axc2 = fig.add_axes((0.81, figy0, 0.1, 0.25), sharey=axc)
        axc3 = axc.twinx()
        axc4 = axc2.twinx()
        axc3.get_shared_y_axes().join(axc3, axc4)
        mpl_br = UTC('2018-05-27').matplotlib_date

        to_plot = np.array(list(zip(mpl, mag, mag**4, mpl)))
        axc.scatter(mpl3[mpl3<=mpl_br], mag3[mpl3<=mpl_br], 4, '0.6', clip_on=False)
        axc2.scatter(mpl3[mpl3>mpl_br], mag3[mpl3>mpl_br], 4, '0.6', clip_on=False)
        axc.scatter(*list(zip(*to_plot[mpl<=mpl_br])), cmap=cmap, norm=norm, linewidths=0.5, edgecolors='k', clip_on=False)
        axc2.scatter(*list(zip(*to_plot[mpl>mpl_br])), cmap=cmap, norm=norm, linewidths=0.5, edgecolors='k', clip_on=False)

        cumM0 = Ml2cumM0(id2, mag2)
        axc3.plot(mpl2, cumM0, 'k')
        cumM0[mpl2<=mpl_br-0.5] = cumM0[mpl2<=mpl_br-0.5][-1]
        axc4.plot(mpl2, cumM0, 'k')
        print(f'sum M0 = {cumM0[-1]:.3e} Nm')

        axc.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
        axc2.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
        axc.set_xlim(UTC('2018-05-09').matplotlib_date, mpl_br-0.5)
        axc2.set_xlim(mpl_br-4, None)

        axc2.tick_params(axis='y', which='both', left=False, labelleft=False)
        axc3.tick_params(axis='y', which='both', right=False, labelright=False)
        plt.setp(axc.get_xticklabels(), rotation=60, ha='right', rotation_mode='anchor')
        plt.setp(axc2.get_xticklabels(), rotation=60, ha='right', rotation_mode='anchor')
        axc4.ticklabel_format(axis='y', useMathText=True)
        axc.set_yticks([1.5, 2, 2.5, 3, 3.5])

        axc.spines['right'].set_visible(False)
        axc2.spines['left'].set_visible(False)
        axc3.spines['right'].set_visible(False)
        axc4.spines['left'].set_visible(False)
        axc.spines['top'].set_visible(False)
        axc2.spines['top'].set_visible(False)
        axc3.spines['top'].set_visible(False)
        axc4.spines['top'].set_visible(False)
        axc4.set_ylim(0, None)
        axc.set_xticks([UTC(f'2018-05-{day}').matplotlib_date for day in (10, 13, 16, 19, 22, 25)])
        axc2.set_xticks([UTC(f'2018-{day}').matplotlib_date for day in ('05-29', '06-08', '06-18')])
        # broken axis
        d = .005
        kw = dict(transform=fig.transFigure, color='k', clip_on=False, lw=1)
        axc.plot((0.8-d, 0.8+d), (figy0-d, figy0+d), **kw)
        axc2.plot((0.81-d, 0.81+d), (figy0-d, figy0+d), **kw)

        axc.set_ylabel(r'WEBNET local magnitude $M_{\rm l}$')
        axc4.set_ylabel('cumulative seismic moment (Nm)')
        akw = dict(xycoords='axes fraction', textcoords='offset points', size='x-large')
        ax3.annotate('a)', (0, 1), (-20, 10), **akw)
        axc.annotate('b)', (0, 1), (-30, 10), **akw)


    if not plot_fm:
        ax4.set_xlabel('depth (km)')
        ax5.set_ylabel('depth (km)')
    if convert_coords:
        if label:
            ax3.set_xlabel('easting (km)')
            ax3.set_ylabel('northing (km)')
            if not plot_fm:
                ax5.set_xlabel('easting (km)')
                ax4.set_ylabel('northing (km)')
    else:
        if label:
            ax3.set_xlabel('longitude')
            ax3.set_ylabel('latitude')
            if not plot_fm:
                ax5.set_xlabel('longitude')
                ax4.set_ylabel('latitude')
    if out:
        plt.savefig(out, dpi=dpi)
    if show:
        plt.show()
    return fig


def plot_events2018(events):
    events1 = zip(*events2lists(events.filter('magnitude > 1.8')))
    all_events = zip(*events2lists(events))
    fig = plot_events_stations_map_depth(events1, convert_coords=LATLON0,
                                         all_events=all_events,
                                         show=False, plot_color_ax=True, colorbar=False)
    fig.axes[0].set_xlim(-1.875, 1.875)
    fig.axes[0].set_ylim(-2.3, 2.7)
    fig.axes[1].set_xlim(5.95, 10.95)
    fig.savefig('../figs/eventmap.pdf', bbox_inches='tight', pad_inches=0.1)


def plot_grond_focal_mechanisms_map(add_inset=True):
    pevents = load_grond()
    events = []
    mts = []
    for evid, ev in pevents.items():
        ref = ev['ref']
        pars = [ref.name, UTC(ref.time), ref.lon, ref.lat, ref.depth/1000,
                ref.magnitude, ev['grond'].moment_tensor]
        events.append(pars)
        mts.append(ev['grond'].moment_tensor)
    fig = plot_events_stations_map_depth(zip(*events), convert_coords=LATLON0,
                                         plot_fm=mts, show=False)
    fig.axes[0].set_xlim(-1.2, 1.2)
    fig.axes[0].set_ylim(-2.3, 2.7)
    fig.axes[0].set_rasterization_zorder(-50)
    if add_inset:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        from plot_focal_mechanisms import plot_grond_focal_mechanisms
        axins = inset_axes(fig.axes[0], width=0.9, height=2.7)
        plot_grond_focal_mechanisms(axins, color_coded=True)
    fig.savefig('../figs/focal_mechanisms_map.pdf', bbox_inches='tight', pad_inches=0.1, dpi=400)


if __name__ == '__main__':
    print('load events')
    allevents = load_full_catalog()
    print('finished loading')
    allevents = obspy.Catalog(sorted(allevents, key=lambda ev: ev.origins[0].time))
    plot_events2018(allevents)
    plot_grond_focal_mechanisms_map()
