#!/usr/bin/python
# -*- coding:utf8 -*-
from __future__ import print_function, division

import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt
import pyraf.iraf as iraf
import glob
import os
import random

path = '/btfidr/data/20130902/007/cam2'
kernel = 'biasEMOFF_A*'
list_of_files = glob.glob(os.path.join(path, kernel))
random.shuffle(list_of_files)
number_of_files = len(list_of_files)

w = width = pyfits.getheader(list_of_files[0])['NAXIS1']
h = height = pyfits.getheader(list_of_files[0])['NAXIS2']
d = depth = number_of_files

index_array = numpy.arange(d)
mean_array = numpy.zeros(d)
std_array = numpy.zeros(d)
cube = numpy.zeros((d,h,w))

cube[0] = pyfits.getdata(list_of_files[0])
for i in range(1, number_of_files):

    cube[i] = pyfits.getdata(list_of_files[i])
    mean = cube[:i].mean(axis=0, dtype=numpy.float64)[800,800]
    std = cube[:i].std(axis=0, dtype=numpy.float64)[800,800]

    mean_array[i] = mean
    std_array[i] = std

    print('{0:02d}\t{1:.2f}\t\t{2:.2f}'.format(i, mean, std))

plt.subplot(2,1,1)
plt.plot(index_array[5:], mean_array[5:], 'ko')
plt.subplot(2,1,2)
plt.plot(index_array[5:], std_array[5:], 'ko')
p = numpy.poly1d(numpy.polyfit(index_array[3:], std_array[3:], deg=2))
plt.show()