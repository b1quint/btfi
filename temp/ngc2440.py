#!/usr/bin/python
"""
This script contains the calcs used to estimate the flux of NGC2440 in data
from Cuesta and Philips 2000.

by Bruno Quint
"""
from __future__ import print_function, division
import aplpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

fits_file = "/home/bquint/Data/n2440.fits"
mask_file = "/home/bquint/Data/NGC2440_cube.mask.fits"

fits_data = pyfits.getdata(fits_file)

fig = aplpy.FITSFigure(fits_file)
fig.show_grayscale()

fig.save("plot.png")
plt.show()
