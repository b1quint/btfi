#/usr/bin/python
# -*- coding: utf8 -*-
import aplpy
import astropy.io.fits as pyfits
import numpy
import numpy.ma as ma
import matplotlib.pyplot as plt

fits_file = "/home/bquint/Data/NGC2440_cube-collapsed.fits"
mask_file = "/home/bquint/Data/NGC2440_cube.mask.fits"
regions = "/home/bquint/Data/NGC2440_Surveys.reg"

mask = pyfits.getdata(mask_file).astype(float)
mask = ma.masked_equal(mask, 0)
mask = numpy.reshape(mask, (mask.shape[0],mask.shape[1],1))
mask = numpy.repeat(mask, 4, axis=2)
mask[:,:,0] *= 1.
mask[:,:,1] *= 0.
mask[:,:,2] *= 0.
mask[:,:,3] *= 0.25

fig = aplpy.FITSFigure(fits_file)
fig.show_grayscale()
plt.imshow(mask)

# fig.show_regions(regions)
fig.tick_labels.set_font(size='small')
fig.save("plot.png")

