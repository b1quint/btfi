#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import print_function, division
import aplpy
import astropy.io.fits as pyfits
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage

__author__ = 'b1quint'
__version__ = '2015.04.06a'

if __name__ == '__main__':

    path = '.'
    path = '/home/bquint/public_html/Data'
    # bname = 'vmap_cubeSAMI_30Dor_Halpha_WCS_med0.fits'  # Background image name
    # cname = 'vmap_cubeVLT_30dor_Halpha_WCS_med0.fits'   # Contour image name
    bname = 'vmap_NGC2440.fits'
    # bname = 'vmap_cube30Dor_SOAR.fits'
    # cname = 'vmap_cube30Dor_VLT.fits'

    bname = os.path.join(path, bname)
    # cname = os.path.join(path, cname)

    fig = aplpy.FITSFigure(bname)

    fig._data = ndimage.median_filter(fig._data, 3)

    vmin = np.median(fig._data) - 0.5 * np.std(fig._data)
    vmax = np.median(fig._data) + 3 * np.std(fig._data)
    # fig._data[fig._data < vmin] = np.nan
    # fig._data[fig._data > vmax] = np.nan
    fig.show_colorscale(cmap='Spectral', vmin=vmin, vmax=vmax)
    fig.set_nan_color('green')
    fig.add_grid()

    # fig.add_colorbar()
    # fig.colorbar.set_axis_label_text('V [km/s]')
    # fig.show_contour(cname, levels=np.linspace(vmin, vmax, 60),
    #                  cmap='Spectral', linewidths=3, returnleevels=True)

    fig.add_scalebar(30./3600.)
    fig.scalebar.set_label('30 "')
    fig.scalebar.set_color('black')

    masked_data = np.ma.masked_array(fig._data, np.isnan(fig._data))
    print(" Background Image: ")
    print("     Median velocity: %.2f km/s" % np.median(masked_data))
    print("     Standard-deviation: %.2f km/s" % np.std(masked_data))
    print("     Minimum: %.2f km/s" % np.min(masked_data))
    print("     Maximum: %.2f km/s" % np.max(masked_data))

    # data = pyfits.getdata(cname)
    # masked_data = np.ma.masked_array(data, np.isnan(data))
    # print(" Contours: ")
    # print("     Median velocity: %.2f m/s" % np.median(masked_data))
    # print("     Standard-deviation: %.2f m/s" % np.std(masked_data))
    # print("     Minimum: %.2f m/s" % np.min(masked_data))
    # print("     Maximum: %.2f m/s" % np.max(masked_data))
    plt.show()