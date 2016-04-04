#!/usr/bin/env python3
# -*- coding: utf8 -*-
"""
    Estimate SAMI/FP Overhead time
    by Bruno Quint
    2016.04.01
"""
import astropy.io.fits as pyfits
import numpy as np
import os

from datetime import datetime, timedelta
from glob import glob

# Input ---
path = '/data/20160402/001'
wildcard = 'cal*.fits'
readout_time = 2.60 # [seconds]

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
exp_time = h['EXPTIME']

t1 = times.min()
t2 = times.max()
dt = t2 - t1
total_time = dt.seconds
total_overhead_time = total_time - exp_time * times.size
total_readout_time = readout_time * times.size
total_fp_overhead_time = total_overhead_time - total_readout_time
average_fp_overhead_time = total_fp_overhead_time / (times.size - 1)

print(" Number of channels: %d" % times.size)
print(" Exposure time per channel: %.2f s" % exp_time)
print("\n Total exposure time: %.2f s" % (times.size * exp_time))
print(" Total time elapsed: %.2f s" % total_time)
print(" Total overhead time: %.2f s" % total_overhead_time)
print("\n Readout time per image: %.2f s" % readout_time)
print(" Total readout time: %.2f s" % total_readout_time)
print(" Total FP overhead time: %.2f s" % total_fp_overhead_time)
print(" Average FP overhead time: %.2f s" % average_fp_overhead_time)
