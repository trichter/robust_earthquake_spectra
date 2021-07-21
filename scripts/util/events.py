# Copyright 2021 Tom Eulenfeld, MIT license

from glob import glob
import grond
import json
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime as UTC
import os.path
from pyrocko.guts import load_all
from pyrocko.model.event import load_events
from pyrocko.moment_tensor import MomentTensor


GRONDRESULTS = '../grond/report/*/*/'
QOPENRESULTS = '../qopen/01_go/results.json'
QOPENEVENTRESULTS = '../qopen/03_mag/results.json'


class MidpointsNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoints=None, clip=False):
        self.midpoints = list(midpoints)
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x = [self.vmin] + self.midpoints + [self.vmax]
        y = np.linspace(0, 1, len(self.midpoints) + 2)
        return np.ma.masked_array(np.interp(value, x, y))


def get_bounds():
    bounds=['2018-05-10', '2018-05-11 03:10:00', '2018-05-12 06:00', '2018-05-21', '2018-05-25', '2018-06-19']
    return [UTC(b) for b in bounds]


def get_bounds2():
    cmap = plt.get_cmap('turbo')
    norm, bounds = get_norm()
    bounds2_mpl= bounds[2:4]
    mpl_values = [np.mean(bounds[:3]), np.mean(bounds[2:4]), np.mean(bounds[3:])]
    colors = [cmap(norm(v)) for v in mpl_values]
    labels = ['before 2018-05-12 06:00', 'period in between', 'after 2018-05-21']
    return bounds2_mpl, colors, labels


def get_norm(extend=False):
    bounds = [b.matplotlib_date for b in get_bounds()]
    if extend:
        bounds[-1]=UTC('2018-05-26').matplotlib_date
    norm = MidpointsNormalize(vmin=bounds[0], vmax=bounds[-1], midpoints=bounds[1:-1])
    return norm, bounds


def event2dict(ev):
    id_ = str(ev.resource_id).split('/')[-1]
    ori = ev.preferred_origin() or ev.origins[0]
    mag = ev.preferred_magnitude() or ev.magnitudes[0]
    try:
        mag = mag.mag
    except:
        mag = None
    picks = None

    return dict(
            id=id_, time=ori.time,
            lat=ori.latitude, lon=ori.longitude, dep=ori.depth / 1000,
            mag=mag, picks=picks)
    pass


def event2list(ev):
    """id, time, lon, lat, dep_km, mag, picks"""
    d = event2dict(ev)
    return (d['id'], d['time'], d['lon'], d['lat'], d['dep'], d['mag'],
            d['picks'])


def events2dicts(events):
    return [event2dict(ev) for ev in events]


def events2lists(events):
    """id, time, lon, lat, dep_km, mag, picks"""
    return [event2list(ev) for ev in events]


def load_grond(load_ensembles=False, load_stats=False):
    globexpr = GRONDRESULTS + 'event.reference.yaml'
    events = {}
    tqdmit = lambda it: it
    if load_ensembles:
        try:
            from tqdm import tqdm as tqdmit
        except ImportError:
            pass
    for fname in tqdmit(glob(globexpr)):
        path = os.path.join(os.path.dirname(fname), '')
        ref = load_events(path + 'event.reference.yaml')[0]
        evid = ref.name[1:]
        mpl = UTC(ref.time).matplotlib_date
        try:
            best = load_events(path + 'event.solution.best.yaml')[0]
        except IndexError:
            print(evid, 'missing in grond results')
            events[evid] = dict(ref=ref, mpl=mpl)
            continue
        mean = load_events(path + 'event.solution.mean.yaml')[0]
        mean_m_dc = mean.moment_tensor.standard_decomposition()[1][2]
        mean.moment_tensor_dc = MomentTensor(m=mean_m_dc)
        events[evid] = dict(ref=ref, mpl=mpl, grond=mean, grond_best=best)
        if load_ensembles:
            ensemble = load_events(path + 'event.solution.ensemble.yaml')
            events[evid]['grond_ensemble'] = ensemble
        if load_stats:
            stats = load_all(filename=path + 'stats.yaml')
            assert len(stats) == 1
            events[evid]['grond_stats'] = stats[0].parameter_stats_list
    events = dict(sorted(events.items(), key=lambda it: it[1]['mpl']))
    return events


def load_qopen(which):
    fname = (QOPENRESULTS if which == 'Q' else
             QOPENEVENTRESULTS if which == 'events' else which)
    with open(fname) as f:
        return json.load(f)


def load_qopen_grond_sds():
    with open('../onsets/M0.json') as f:
        M0s_onset = json.load(f)
    events = load_grond()
    events2 = {}
    q = load_qopen('events')
    for evid, ev in events.items():
        if evid not in q['events'] or 'Mcat' not in q['events'][evid]:
            print(evid, 'missing in qopen results')
            continue
        ev['qopen'] = q['events'][evid]
        ev['M0_onset'] = np.median(M0s_onset[evid])
        events2[evid] = ev
    return events2
