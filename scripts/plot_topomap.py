# Copyright 2021 Tom Eulenfeld, MIT license

from cartopy.feature import NaturalEarthFeature as NEF
import glob
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
import obspy
import shapely.geometry as sgeom

from util.topomap import EA_EURO, GEO, add_scale, de_border, add_ticklabels, plot_elevation
from util.events import events2lists
from prepare_data import load_full_catalog


def load_webnet(year='*', plot=False, stations=False):
    try:
        import pandas as pd
    except ImportError:
        print('pandas is needed to plot WEBNET catalog')
        eqs = None
    else:
        names = ['time', 'lon', 'lat', 'dep', 'mag']
        kwargs = dict(sep=';', skipinitialspace=True, skiprows=3, parse_dates=[0],
                      names=names)
        path ='../data/webnet/*.txt'
        frames = [pd.read_csv(fname, **kwargs) for fname in glob.glob(path)]
        if len(frames) == 0:
            print('You can obtain the WEBNET catalog at ig.cas.cz. '
                  'Put the txt files in new webnet directory inside data.')
            eqs = None
        else:
            eqs = pd.concat(frames, ignore_index=True)
    if stations:
        path = '../data/station_coordinates.txt'
        sta = pd.read_csv(path, sep='\s+', usecols=(0, 1, 2))
    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.scatter(eqs.lon.values, eqs.lat.values, 1, eqs.time.values)
        ax2.scatter(eqs.time.values, eqs.mag.values, 1, eqs.time.values)
        if stations:
            ax1.scatter(sta.lon.values, sta.lat.values, 100, marker='v',
                        color='none', edgecolors='k')
            for idx, s in sta.iterrows():
                ax1.annotate(s.station, (s.lon, s.lat), (5, 5),
                             textcoords='offset points')
        plt.show()
    if stations:
        return eqs, sta
    return eqs


def plot_topomap(events=None):
    mlf = [(12.410, 50.261), (12.485, 50.183), (12.523, 50.127),
           (12.517, 50.125), (12.538, 50.131), (12.534, 50.130),
           (12.543, 50.109), (12.547, 50.083), (12.547, 50.074),
           (12.545, 50.066), (12.546, 50.057), (12.576, 50.034),
           (12.594, 50.016), (12.632, 50.004), (12.665, 49.980)]
    er1 = [(12.491, 50.167), (12.522, 50.180), (12.598, 50.205), (12.688, 50.243)]
    er2 = [(12.554, 50.118), (12.585, 50.128), (12.623, 50.143), (12.679, 50.161)]
    wfz = [(12.507, 50.217), (12.562, 50.180), (12.597, 50.149), (12.625, 50.132)]
    tfz1 = [(12.249, 50.174), (12.213, 50.196), (12.200, 50.210), (12.188, 50.221), (12.165, 50.243)]
    tfz2 = [(12.318, 50.095), (12.358, 50.079), (12.364, 50.074), (12.421, 50.011), (12.439, 49.997)]
    ufz = [(12.626, 50.235), (12.651, 50.245), (12.667, 50.248), (12.676, 50.252)]
    novy_kostel = (12.4210, 50.2203)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=EA_EURO)
    extent = [12.19, 12.63, 50.05, 50.35]
    ax.set_extent(extent, crs=GEO)

    box = ax.get_position().bounds
    subax = fig.add_axes([box[0]-0.05, box[1]+box[3]-0.3, 0.28, 0.28], projection=EA_EURO)
    subax.set_extent([8, 16, 47, 55], GEO)
    subax.add_feature(NEF('physical', 'land', '10m'), facecolor='0.7', alpha=0.5, rasterized=True)
    subax.add_feature(NEF('physical', 'coastline', '10m'), facecolor='none', edgecolor='k', linewidth=0.5, rasterized=True)
    subax.add_feature(NEF('cultural', 'admin_0_boundary_lines_land', '10m'), facecolor='none', edgecolor='k', linewidth=0.5, rasterized=True)
    subax.add_geometries([sgeom.box(extent[0], extent[2], extent[1], extent[3])], GEO,
                          facecolor='none', edgecolor='k', linewidth=1, alpha=0.8)
    subax.annotate('Germany', (9.5, 50), None, GEO._as_mpl_transform(subax), size='small', color='0.5', rotation=70)
    subax.annotate('Czechia', (12.6, 49.5), None, GEO._as_mpl_transform(subax), size='x-small', color='0.5', rotation=10)
    lonticks = [12.2, 12.3, 12.4, 12.5, 12.6]
    latticks = [50.1, 50.2, 50.3]
    add_ticklabels(ax, lonticks, latticks)
    ax.tick_params(axis='both', which='major', labelsize=8)
    plot_elevation(ax, shading=False, cmap=LinearSegmentedColormap.from_list('littlegray', ['white', '0.5']),
                   azimuth=315, altitude=60, rasterized=True)
    add_scale(ax, 5, (12.47, 50.05), color='0.2')
    de_border(ax, edgecolor='0.6', rasterized=True)
    eqs, sta = load_webnet(stations=True)
    if eqs is not None:
        ax.scatter(eqs.lon.values, eqs.lat.values, 4, '#9bd7ff', alpha=0.4, marker='.', transform=GEO, rasterized=True)
    if events is not None:
        _, _, lon, lat, _, mag, *_ = zip(*events2lists(events))
        ax.scatter(lon, lat, 4, 'C0', marker='o', transform=GEO)
    ax.plot(*list(zip(*mlf)), color='0.5', transform=GEO)
    for lineament in (er1, er2, wfz, tfz1, tfz2, ufz):
        ax.plot(*list(zip(*lineament)), color='0.5', transform=GEO, lw=1)
    geotrans = GEO._as_mpl_transform(ax)
    ax.annotate('MLF', (12.505, 50.145), None, geotrans, size='x-small', zorder=10, rotation=295)
    ax.annotate('Eger Rift', (12.55, 50.15), None, geotrans, color='0.2', zorder=10, rotation=30)
    ax.annotate('Cheb\n   Basin', (12.4, 50.11), None, geotrans, color='0.2', zorder=10, rotation=-10)
    used_stations = 'NKC LBC VAC KVC STC POC SKC KRC ZHC'.split()
    sta = sta[sta.station.isin(used_stations)]
    ax.scatter(sta.lon.values, sta.lat.values, 100, marker='^', color='none', edgecolors='k', transform=GEO, zorder=10)
    ax.scatter(*novy_kostel, 20, marker='o', color='none', edgecolors='k', transform=GEO, zorder=10)
    ax.annotate('Nový\nKostel', novy_kostel, (-5, 0), geotrans, 'offset points', size='x-small', ha='right', va='top', zorder=10)
    for idx, s in sta.iterrows():
        xy = (2, 2) if s.station not in ('KOPD', 'KAC') else (-10, 5)
        ax.annotate(s.station, (s.lon, s.lat), xy, geotrans, 'offset points', size='x-small', zorder=10)
    x0, y0 = EA_EURO.transform_point(LATLON0[1], LATLON0[0], GEO)
    ax.add_geometries([sgeom.box(x0-1875, y0-2300, x0+1875, y0+2700)], EA_EURO,
                       facecolor='none', edgecolor='C1', linewidth=2, alpha=0.8, zorder=11)
    fig.savefig('../figs/topomap.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
    return sta


LATLON0 = (50.25, 12.45)
MAG = 1.5

if __name__ == '__main__':
    print('load events')
    allevents = load_full_catalog()
    print('finished loading')
    allevents = obspy.Catalog(sorted(allevents, key=lambda ev: ev.origins[0].time))
    events = allevents.filter(f'magnitude > {MAG}')
    sta = plot_topomap(events=allevents)
