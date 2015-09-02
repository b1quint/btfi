#/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import division, print_function

import numpy
import matplotlib.pyplot as plt
import os
import scipy.fftpack as fftpack

path = '/home/bquint/Dropbox/BTFI (1)/Fabry-Perot/data'
filename = 'cs100uofh_10hz2_S1C.csv'

t, v = numpy.genfromtxt(os.path.join(path, filename), delimiter=',', dtype=str, unpack=True)
t, v = t[2:].astype(float), v[2:].astype(float)

fv = fftpack.fft(v)
# sfv = fftpack.fftshift(fv)
ft = fftpack.fftfreq(fv.size, d=(t[1:] - t[:-1]).mean())

plt.subplot(2,1,1)
plt.plot(t, v)

plt.subplot(2,1,2)
plt.plot(ft[ft > 0], numpy.abs(fv[ft > 0]))
plt.xlim(xmin=0, xmax=500)
plt.grid()
plt.show()