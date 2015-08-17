#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    Measure instrument readout time

    This script estimates the readout time based on the header's information.
    It basically creates a table that gets the time that each file was created
    and compares with the exposure time.

    The input arguments are the names of the files so try not to mix files that
    have different exposure times.
"""

from __future__ import print_function, division

import astropy.io.fits as pyfits
import glob
import numpy
import os
from datetime import datetime

if __name__ == '__main__':

    path = "/btfidr/data/20130824/019/cam2"
    file_pattern = "darkD*.fits"
    files = glob.glob(os.path.join(path, file_pattern))
    files.sort()

    exp_time = []
    shutter_time = []
    for filename in files:
        header = pyfits.getheader(filename)
        exp_time.append(float(header['exptime']))
        shutter_time.append(datetime.strptime(header['utshut'], "%Y-%m-%dT%H:%M:%S.%f"))

    exp_time = numpy.array(exp_time)
    exp_time = numpy.average(exp_time)

    shutter_time = numpy.array(shutter_time)
    mean_exp_time = shutter_time[1:] - shutter_time[:-1]
    mean_exp_time = mean_exp_time.astype('timedelta64[ms]').astype('int')
    mean_exp_time = numpy.average(mean_exp_time) / 1000.

    print('Total read time: {0:.2f}'.format(mean_exp_time - exp_time))

