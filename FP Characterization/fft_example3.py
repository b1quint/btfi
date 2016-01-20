#!/usr/bin/python
# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
import numpy

from scipy import fftpack
from scipy import signal

dt = 0.001
t = numpy.arange(0, 1, dt)
x = numpy.zeros_like(t)
y = numpy.abs(signal.sawtooth(2 * numpy.pi * 100 * t)) - 0.5
y = numpy.where(y > 0, y, 0)
x += y
x += signal.square(2 * numpy.pi * 10 * t)

plt.subplot(2,1,1)
plt.plot(t, x)
plt.grid()

fx = fftpack.fft(x)
ft = fftpack.fftfreq(fx.size, d=dt)

plt.subplot(2,1,2)
plt.plot(ft,numpy.abs(fx))
plt.grid()
plt.show()
