#!/usr/bin/env python3
# -*- coding: utf8 -*-
"""
    Estimate SAMI's readout time
    by Bruno Quint
    2016.04.01
"""
import astropy.io.fits as pyfits
import numpy as np
import os

from datetime import datetime
from glob import glob

# Input ---
path = '/data/20160402/001'
wildcard = '*.fits'

# Main thread ---
files = glob(os.path.join(path, wildcard))
files = np.array(files)
times = []
for f in files:
    h = pyfits.getheader(f)
    t = datetime.strptime("%s %s" % (h['DATE-OBS'], h['TIME-OBS']),
            "%Y-%m-%d %H:%M:%S.%f")
    times.append(t)

times = np.array(times)
t1 = times.min()
t2 = times.max()
dt = t2 - t1
readout_time = dt.seconds / times.size
print(" Estimated read-out time: %.2f seconds" % readout_time)
