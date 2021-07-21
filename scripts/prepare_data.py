# Copyright 2021 Tom Eulenfeld, MIT license

import os
import shutil
from glob import glob
import pickle

from obspy import read_events, read_inventory
from pyrocko.model.event import dump_events
from pyrocko.io.quakeml import QuakeML
from pyrocko.gui.marker import PhaseMarker


DATA = '../data/'
TMP = '/tmp/'
PHA_Q1 = DATA + 'catalog_2018swarm_quality1.pha'
PHA_ALL = DATA + 'catalog_2018swarm.pha'
QUAKEMLOUT = TMP + 'catalog_2018_qu1.xml'
PYROCKOEVENTOUT = DATA + 'pyrocko_events/'
INV = DATA + 'webnet_9stations.xml'
WFPATH1 = DATA + 'waveforms/{evid}_*.mseed'
WFPATH2 = DATA + 'waveforms_event/x{evid}/'
PKL_EVENTS = TMP + 'catalog_2018.pkl'


def load_full_catalog():
    if not os.path.exists(PKL_EVENTS):
        events = read_events(PHA_ALL)
        with open(PKL_EVENTS, 'wb') as f:
            pickle.dump(events, f, pickle.HIGHEST_PROTOCOL)
    with open(PKL_EVENTS, 'rb') as f:
        return pickle.load(f)


def _filter_pha():
    qcfile = DATA + '2018_events_qc_mag1.9.txt'
    with open(qcfile) as f:
        text = f.read()
    cat = load_full_catalog()
    print('All events', cat)
    qc = {line.split()[0]: line.split(maxsplit=2)[2] for line in text.splitlines() if line.startswith('2018')}
    events = []
    for event in cat.events:
        id_ = str(event.resource_id).split('/')[-1]
        try:
            fields = qc[id_].split()
        except KeyError:
            continue
        if 'Q1' in fields:
            events.append(event)
    cat.events = events
    print('Only Q1 events', cat)
    cat.write(PHA_Q1, 'HYPODDPHA')


def copywf(evid):
    path1 = WFPATH1.format(evid=evid)
    path2 = WFPATH2.format(evid=evid)
    if not os.path.exists(path2):
        os.makedirs(path2)
        for fname in glob(path1):
            print(f'copy file {fname}')
            shutil.copy(fname, path2)


def convert_to_pyrocko():
    if not os.path.exists(PYROCKOEVENTOUT):
        os.mkdir(PYROCKOEVENTOUT)
    inv = read_inventory(INV)
    obs_events = read_events(PHA_Q1, inventory=inv)
    obs_events.write(QUAKEMLOUT, 'QUAKEML')
    qml = QuakeML.load_xml(filename=QUAKEMLOUT)
    events = qml.get_pyrocko_events()
    picks = qml.get_pyrocko_phase_markers()
    for event in events:
        event.name = 'x' + event.name.split('/')[-1]
        copywf(event.name[1:])
        path = PYROCKOEVENTOUT + event.name
        dump_events([event], path + '_event.pf')
        eventpicks = [p for p in picks if p.get_event().name.split('/')[-1] == event.name[1:]]
        for p in eventpicks:
            p.set_event(event)
        PhaseMarker.save_markers(eventpicks, path + '_picks.txt')


# the following two functions are not used, because we only use WEBNET stations
from obspy.clients.fdsn import Client

INV2 = DATA + 'webnet_sx.xml'
WFPATH2_SX = DATA + 'waveforms_event/x{evid}/{evid}_SX_{station}_{channel}.mseed'


def download_SX_stations():
    client = Client('BGR')
    inv2 = client.get_stations(
        starttime='2018-01-01', endtime='2019-01-01',
        network='SX', channel='HH?', level='response',
        latitude=50.25, longitude=12.45, maxradius=18/111)
    inv = read_inventory(INV) + inv2
    inv.write(INV2, 'STATIONXML')


def download_SX_waveforms():
    inv = read_inventory(INV2).select(network='SX')
    client = Client('BGR')
    obs_events = read_events(PHA_Q1)
    for event in obs_events:
        evid = str(event.resource_id).split('/')[-1]
        otime = event.origins[0].time
        print('event', evid)
        for sta in inv[0]:
            print('  download waveforms for station', sta.code)
            stream = client.get_waveforms('SX', sta.code, '', 'HH?',
                                          otime-10, otime+20)
            for tr in stream:
                path2 = WFPATH2_SX.format(evid=evid, station=sta.code,
                                          channel=tr.stats.channel)
                tr.write(path2, 'MSEED')


if __name__ == '__main__':
    _filter_pha()
    convert_to_pyrocko()
