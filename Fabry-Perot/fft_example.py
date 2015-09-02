#/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import division, print_function

import numpy
import matplotlib.pyplot as plt
import os
import scipy.fftpack as fftpack

def add_frequency(t, f):
    return numpy.sin(f * t * (2 * numpy.pi))

dt = 0.001
t = numpy.arange(0, 10, dt)
x = numpy.zeros_like(t)
x += add_frequency(t, 1)
x += add_frequency(t, 3)

plt.subplot(2,1,1)
plt.plot(t, x)
plt.grid()

fx = fftpack.fft(x)
ft = fftpack.fftfreq(fx.size, d=dt)

plt.subplot(2,1,2)
plt.plot(ft,numpy.abs(fx))
plt.grid()
plt.show()