# Copyright 2021 Tom Eulenfeld, MIT license

import numpy as np

def convert_coords2km(coords, latlon0=None):
    import utm
    x, y = zip(*[utm.from_latlon(lat1, lon1)[:2]
                 for lat1, lon1, *_ in coords])
    if latlon0 is None:
        x0 = np.mean(x)
        y0 = np.mean(y)
    else:
        x0, y0 = utm.from_latlon(*latlon0)[:2]
    x = (np.array(x) - x0) / 1000
    y = (np.array(y) - y0) / 1000
    if len(coords[0]) == 3:
        return list(zip(x, y, [c[2] for c in coords]))
    else:
        return list(zip(x, y))
