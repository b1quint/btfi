#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function, division

import astropy.io.fits as pyfits
import glob
import numpy
import os
from datetime import datetime

if __name__ == '__main__':

    readout_time = 1.15  # seconds
    step_time = 2.86  # seconds

    path = "/btfidr/data/20130824/015/cam2"
    file_pattern = "scanB*.fits"
    files = glob.glob(os.path.join(path, file_pattern))
    files.sort()

    exp_time = []
    shutter_time = []
    for filename in files:
        header = pyfits.getheader(filename)
        exp_time.append(float(header['exptime']))
        shutter_time.append(datetime.strptime(header['utshut'], "%Y-%m-%dT%H:%M:%S.%f"))
        print(filename, header['exptime'], header['utshut'])

    exp_time = numpy.array(exp_time)
    exp_time = numpy.average(exp_time)

    shutter_time = numpy.array(shutter_time)
    mean_overhead_time = shutter_time[1:] - shutter_time[:-1]
    mean_overhead_time = mean_overhead_time.astype('timedelta64[ms]').astype('int')
    mean_overhead_time = mean_overhead_time / 1000. - exp_time - readout_time
    mean_overhead_time = numpy.average(mean_overhead_time[mean_overhead_time > 0.1])

    print('Total overhead time: {0:.2f}'.format(mean_overhead_time))