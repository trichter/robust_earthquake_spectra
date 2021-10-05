# Copyright 2021 Tom Eulenfeld, MIT license

import json
import matplotlib.pyplot as plt
from mtspec import mtspec
import numpy as np
from obspy import read, read_events, read_inventory
from obspy.core.event.source import farfield
from obspy.geodetics import gps2dist_azimuth

from util.events import load_grond, load_qopen


OUT = '../figs/'
DATA = '../data/'
STATIONS = 'KRC KVC LBC NKC POC SKC STC VAC ZHC'.split()


def dist(inv, stream, event):
    coords = inv.get_coordinates(stream[0].id)
    ori = event.origins[0]
    coord_args = (coords['latitude'], coords['longitude'], ori.latitude, ori.longitude)
    d, az1, az2 = gps2dist_azimuth(*coord_args)
    inc = np.arctan2(d, ori.depth)
    # event to station coordinates in North-East-Down convention
    x = -np.cos(np.deg2rad(az1)) * d
    y = -np.sin(np.deg2rad(az1)) * d
    z = -ori.depth
    d = (d ** 2 + ori.depth ** 2) ** 0.5
    print(stream[0].id, *np.round(coord_args, 2), x, y, z)
    return d, np.deg2rad(az1), inc, (x, y, z)


def get_picktime(picks, station, phase):
    picks = [p for p in picks if p.waveform_id.station_code == station and
             p.phase_hint == phase]
    assert len(picks) == 1
    return picks[0].time


def load_displacement_data(evid, sta, phase='P'):
    stream = read(DATA + f'waveforms_event/x{evid}/{evid}_{sta}*.mseed')
    stream = stream.select(station=sta)
    stream.detrend('demean')
    stream.remove_response(INV, 'DISP')
    stream.detrend('demean')
    event = [e for e in EVENTS if str(e.resource_id).split('/')[-1] == evid][0]
    d, az, inc, xyz = dist(INV, stream, event)
    t = get_picktime(event.picks, sta, phase)
    if phase == 'P':
        tr = stream.select(component='Z')[0]
    else:
        stream.trim(t+TW[0]-1, t+TW[1]+1)
        stream.rotate('NE->RT', az)
        tr = stream.select(component='T')[0]
    tr.trim(t+TW[0], t+TW[1])
    return tr, d, xyz


def ampl2sds(spec, dist, rp=1, fs=2, rho=2700, phase='P'):
    # Hanks Wyss 1972
    v = 5780 if phase == 'P' else 3400
    return spec / fs / rp * 4 * np.pi * rho * dist * v ** 3


def damping(dist, q, at_freq=None, phase='S'):
    # simply assume Qp=Qs
    # v = 5780 if phase == 'P' else 3400
    v = 3400
    geff = np.array(q['b']) / v + np.array(q['g0'])
    if at_freq is not None:
        geff = np.interp(at_freq, q['freq'], geff)
    damp = np.exp(-dist / 2 * geff)
    print(f'low freq damping for {phase} wave: {damp[0]:.2f}')
    return damp


def radpat(mt, xyz, phase, norm=1):
    if xyz is None:
        theta = np.arange(0, 181, 1)
        phi = np.arange(0, 360, 1)
        theta, phi = np.meshgrid(theta, phi)
        tp = np.vstack([np.ravel(theta), np.ravel(phi)])
        disp2 = farfield(mt.m6(), tp, phase)
        disp2 = np.sum(disp2**2, axis=0) ** 0.5
        disp = np.max(disp2)
        return disp
    xyz = np.array(xyz)[:, np.newaxis]
    disp = farfield(mt.m6(), xyz, phase)
    disp = (np.sum(disp**2)) ** 0.5
    rp = disp / norm
    print(f'radiation pattern: {disp:.2e}/{norm:.2e}={rp:.3}   moment: {mt.moment:.2e}')
    return rp


def calc_onset_spectra_and_compare(evid, compare_without_radpat=False,
                                   correct_damping=True, ext='pdf',
                                   label=None, savefig=True):
    q = load_qopen('events')
    q2 = load_qopen('Q')
    try:
        mt = load_grond()[evid]['grond'].moment_tensor
    except KeyError:
        print(evid, 'missing in grond results')
        return
    qfreq = q['freq']
    try:
        qsds = q['events'][evid]['sds']
        qM0 = q['events'][evid]['M0']
        qfc = q['events'][evid]['fc']
    except KeyError:
        print(evid, 'missing in qopen results')
        return
    fig = plt.figure()
    ax1 = fig.add_subplot(121 if compare_without_radpat else 111)
    onset_sds = []
    onset_M0s = []
    for phase in 'PS':
        dispm = radpat(mt, None, phase)
        for i, sta in enumerate(STATIONS):
            print(sta, phase)
            try:
                tr, d, xyz = load_displacement_data(evid, sta, phase)
            except Exception as ex:
                print(ex)
                continue
            print(f'distance {d:.0f}m')
            ampl_spec, freq = mtspec(
                    tr.data, tr.stats.delta, 2)
            ampl_spec=ampl_spec[freq>0]
            freq=freq[freq>0]
            p = '-' if phase == 'P' else '--'
            rp = radpat(mt, xyz, phase, norm=dispm)
            spec = ampl2sds(ampl_spec**0.5, d, rp=rp, phase=phase)
            if correct_damping:
                spec = spec / damping(d, q2, freq, phase=phase)
            ax1.plot(freq, spec, p, color=str((i+5)/ 16), label=f'{sta} {phase}')
            onset_sds.append([phase, sta, np.round(freq, 2).tolist(),
                              np.round(spec).tolist()])
            M0 = np.round(np.mean(spec[freq<3.5]))
            onset_M0s.append(M0)
    for _ in range(4):
        ax1.plot([], [], 'wx', label=' ')  # ghost legend entry
    ax1.plot(qfreq, qsds, color='C1', lw=2, label='Qopen')
    ax1.axhline(qM0, 0.6, 0.89, color='C1', lw=2, ls='--', label='Qopen $M_0$', zorder=12)
    fc_handle = ax1.axvline(qfc, 0.5, color='C1', lw=2, ls='--', label=r'Qopen $f_{\rm c}$', alpha=0.5, zorder=-2)
    ax1.axhline(mt.moment, 0.6, 0.95,  color='C0', lw=2, ls='--', label='Grond $M_0$', zorder=11)
    ax1.axhline(np.median(onset_M0s), 0.6, color='C6', lw=2, ls='--', label='onset $M_0$', zorder=10)
    ax1.annotate(evid, (0.98, 0.98), xycoords='axes fraction',
                 ha='right', va='top')
    if label:
        ax1.annotate(label, (0.02, 0.9), xycoords='figure fraction', size='large')
    if compare_without_radpat:
        ax2 = fig.add_subplot(122, sharex=ax1, sharey=ax1)
        for phase in 'PS':
            for i, sta in enumerate(STATIONS):
                try:
                    tr, d, xyz = load_displacement_data(evid, sta, phase)
                except Exception as ex:
                    print(ex)
                    continue
                print(f'distance {d:.0f}m')
                spec, freq = mtspec(
                        tr.data, tr.stats.delta, 2)
                spec=spec[freq>0]
                freq=freq[freq>0]
                p = '-' if phase == 'P' else '--'
                amp = ampl2sds(spec**0.5, d, rp=0.5, phase=phase)
                ax2.plot(freq, amp, p, color=str((i+5)/ 16), label=f'{sta} {phase}')
    ax1.set_xlim(qfreq[0], qfreq[-1])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks((1, 10, 100))
    ax1.set_xticklabels(('1', '10', '100'))
    ax1.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'source displacement spectrum $\omega M$ (Nm)')

    # https://stackoverflow.com/a/53126116
    from matplotlib.legend_handler import HandlerLine2D
    def update_prop(handle, orig):
        handle.update_from(orig)
        x,y = handle.get_data()
        handle.set_data([np.mean(x)]*2, [-1.5*y[0], 3.5*y[0]])
    ax1.legend(ncol=3, frameon=False,
               handler_map={fc_handle:HandlerLine2D(update_func=update_prop)})
    if savefig:
        fig.savefig(OUT + f'source_spectrum_{evid}_tw{TW[0]}s_{TW[1]}s.{ext}',
                    bbox_inches='tight', pad_inches=0.1)
    return onset_M0s, onset_sds


def plot_sds(allsds=True, correct_damping=True, savefig=True):
    global INV, EVENTS, TW
    INV = read_inventory(DATA +'webnet_9stations.xml')
    EVENTS = read_events(DATA + 'catalog_2018swarm_quality1.pha', inventory=INV)
    TW = (-0.1, 0.5)

    if allsds:
        M0 = {}
        sds = {}
        for event in EVENTS:
            evid = str(event.resource_id).split('/')[-1]
            res = calc_onset_spectra_and_compare(evid, ext='png', correct_damping=correct_damping, savefig=savefig)
            if res is not None:
                M0[evid], sds[evid] = res
        adds =  '' if correct_damping else '_uncorrected'
        with open(f'../onsets/M0{adds}.json', 'w') as f:
            json.dump(M0, f)
        with open(f'../onsets/sds{adds}.json', 'w') as f:
            json.dump(sds, f)
    else:
        calc_onset_spectra_and_compare('20186784', label='a)', correct_damping=correct_damping, savefig=savefig)
        calc_onset_spectra_and_compare('201856087', label='b)', correct_damping=correct_damping, savefig=savefig)


def _v0_estimate():
    q1 = load_qopen('Q')
    q2 = load_qopen('events')
    v0s = [ev['v0'] for ev in q1['events'].values()]
    print('v0', np.mean(v0s))
    v0s = [ev['v0'] for ev in q2['events'].values()]
    print('v0', np.mean(v0s))


if __name__ == '__main__':
    plot_sds()
    plot_sds(correct_damping=False, savefig=False)
    plot_sds(allsds=False)
