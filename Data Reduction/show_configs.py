#!/usr/bin/python
# -*- coding: utf8 -*-
"""Show several frames of a single data-cube using only one colorbar."""

from __future__ import print_function, division

__author__ = 'Bruno Quint'
__date__ = '2015.03.09'

import numpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.io.fits as pyfits


def main():
    """Main method."""

    # Setup
    filename = '/data/BTFI/20140402/CALIBRATION/cubeD_Neon_6600.5A-16.3A_2014Apr02.fits'
    nframes = 5

    # Preparing figure
    cm2inches = lambda x, y: (x / 2.54, y / 2.54)
    my_fig = plt.figure(figsize=cm2inches(20,5), dpi=96)
    gs = gridspec.GridSpec(2,5, height_ratios=(20,1))

    # Loading data
    data = pyfits.getdata(filename)
    header = pyfits.getheader(filename)
    depth = data.shape[0]

    # Plotting
    vmin = data.mean() - data.std()
    vmax = data.mean() + 2 * data.std()
    config = {'origin': 'lower',
              'interpolation': 'nearest',
              'vmin': vmin,
              'vmax': vmax,
              'cmap': 'gray'}

    axes = []
    for i in range(nframes):
        print(i * depth // nframes)
        ax = my_fig.add_subplot(gs[:-1,i])
        im = ax.imshow(data[i * 5], **config)
        ax.set_xlabel("z = %d" % (i * header['C3_3']))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid()
        axes.append(ax)

    # Display image
    cbar = plt.colorbar(im, cax=plt.subplot(gs[1,:]), orientation='horizontal')
    plt.tight_layout()
    plt.show()

    return


if __name__ == '__main__':
    main()