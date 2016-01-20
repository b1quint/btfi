#/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import division, print_function

import numpy
import matplotlib.pyplot as plt
import os

from scipy import signal
from scipy import fftpack

# t = numpy.linspace(0, 1, 500, endpoint=False)
# plt.plot(t, signal.square(2 * np.pi * 5 * t))
# plt.ylim(-2, 2)

def add_frequency(t, f):
    return numpy.sin(f * t * (2 * numpy.pi))

dt = 0.001
t = numpy.arange(0, 0.5, dt)
x = numpy.zeros_like(t)
x += signal.square(2 * numpy.pi * 10 * t)
x += numpy.abs(numpy.sin(2 * numpy.pi * 50 * t))
# x += signal.square(2 * numpy.pi * 200 * t)
# x += signal.square(2 * numpy.pi * 200 * t)
# x += numpy.sin(2 * numpy.pi * 100 * t)
# x += numpy.sin(2 * numpy.pi * 200 * t)

plt.subplot(2,1,1)
plt.plot(t, x, lw=2)
plt.grid()

fx = fftpack.fft(x)
ft = fftpack.fftfreq(fx.size, d=dt)

plt.subplot(2,1,2)
plt.plot(ft,numpy.abs(fx))
plt.xlim(xmin=0, xmax=500)
plt.grid()
plt.show()