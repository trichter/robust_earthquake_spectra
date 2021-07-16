from obspy import read
DATA = '../data/waveforms_event/x{}/{}_{}_{}.mseed'

def get_waveforms(event=None, station=None, channel=None, starttime=None, endtime=None, **kw):
    id_ = event.resource_id.id.split('/')[-1]
    channel = 'H' + channel[1:]  # file names do not fit the channel code
    stream = read(DATA.format(id_, id_, station, channel), 'MSEED')
    stream.trim(starttime, endtime)
    return stream