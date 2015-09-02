#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function, division

import astropy.io.fits as pyfits
import glob
import numpy
import os
from datetime import datetime, timedelta

if __name__ == '__main__':

    readout_time = 1.15  # seconds

    # path = "/btfidr/data/20130824/012/cam2"
    # file_pattern = "scanD*.fits"
    path = "/btfidr/data/20130824/017/cam2"
    file_pattern = "scanG*.fits"

    files = glob.glob(os.path.join(path, file_pattern))
    files.sort()

    time_start = 0.0
    time_stop = 0.0
    exp_time = []
    shutter_time = []

    for filename in files:
        header = pyfits.getheader(filename)
        exp_time.append(float(header['exptime']))
        shutter_time.append(datetime.strptime(header['utshut'], "%Y-%m-%dT%H:%M:%S.%f"))

    exp_time = numpy.array(exp_time)
    exp_time = numpy.average(exp_time)

    shutter_time = numpy.array(shutter_time)
    mean_overhead_time = shutter_time[1:] - shutter_time[:-1]
    mean_overhead_time = mean_overhead_time.astype('timedelta64[ms]').astype('int')
    mean_overhead_time = numpy.average(mean_overhead_time) / 1000. - readout_time - exp_time

    shutter_time = shutter_time.astype('datetime64[ms]').astype('int')
    shutter_time /= 1000.
    tot_elapsed_time = shutter_time.max() - shutter_time.min()
    tot_elapsed_time = timedelta(seconds=tot_elapsed_time)

    nimages = shutter_time.size
    nsteps = shutter_time.size / float(header['IMGPPOS']) - 1

    print('\nTotal ellapsed time:', tot_elapsed_time)
    print('Exposure time per image: {0:.2f}'.format(exp_time))
    print('Readout time: {0:.2f}'.format(readout_time))
    print('\nNumber of images taken: {0:d}'.format(nimages))
    print('Number of steps given: {0:.0f}'.format(nsteps))
    print('Number of images per position: {0:s}'.format(header['IMGPPOS']))
    print('\nTotal exposure time:', timedelta(seconds=nimages * exp_time))
    print('Total time lost due readout:', timedelta(seconds=nimages * readout_time))
    print('Total time lost in steps:', timedelta(seconds=nsteps * mean_overhead_time))
    print('Checking total time:', timedelta(seconds=(nsteps * mean_overhead_time + nimages * readout_time + nimages * exp_time)))