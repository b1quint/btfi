#!/usr/bin/python
# -*- coding:utf8 -*-

from __future__ import print_function, division

import numpy
import matplotlib.pyplot as plt
import pyraf.iraf as iraf
import glob
import os

iraf.images()
iraf.immatch()

option = 'average'  # 'average'|'median'

path = '/btfidr/data/20130902/007/cam2'
kernel = 'biasEMOFF_A*'
list_of_files = glob.glob(os.path.join(path, kernel))
list_of_files.sort()
number_of_files = len(list_of_files)

output_pattern = '%03d-mean.fits' if option is 'average' else '%03d-midpt.fits'
logfile = 'imstat.log.mean'if option is 'average' else 'imstat.log.midpt'

for i in range(2, number_of_files):
    print(i)
    filename = "list_%03d.temp.dat" % i

    if os.path.exists(filename):
        os.remove(filename)

    if os.path.exists(output_pattern%i):
        continue

    foo = open(filename, 'w')
    foo.write("\n".join(list_of_files[:i]))
    foo.write("\n")
    foo.close()

    iraf.imcombine(input='@%s'%filename, output=output_pattern%i,
                   Stdout='imcombine.stdout', combine=option)
    os.remove(filename)

if option == 'average':
    iraf.imstat(images='*-mean.fits', fields='image,mean,midpt,mode,stddev',
                Stdout=logfile)
else:
    iraf.imstat(images='*-midpt.fits', fields='image,mean,midpt,mode,stddev',
                Stdout=logfile)

if os.path.exists('imcombine.stdout'):
    os.remove("imcombine.stdout")

# mean_file = 'imstat.log.mean'
# mid_file = 'imstat.log.midpt'
#
# std_mean = numpy.loadtxt(mean_file, usecols=(4,), unpack=True)
# std_mid = numpy.loadtxt(mid_file, usecols=(4,), unpack=True)
# x = numpy.arange(std_mean.size)+2
#
# plt.plot(x, std_mean, 'kx')
# plt.plot(x, std_mid, 'ko')
# plt.plot()
# plt.grid()
# plt.xlabel("Number of images combined")
# plt.ylabel("Standard Deviation [counts]")
# plt.show()