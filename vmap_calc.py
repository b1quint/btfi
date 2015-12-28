#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import print_function, division
import astropy.io.fits as pyfits
import os
import numpy as np

from scipy.optimize import leastsq
from scipy import stats
from numpy import ma

__author__ = 'b1quint'
__version__ = '20150406a'

class Main():

    def __init__(self, path, filename):

        self.path = path
        self.filename = filename
        return

    def run(self):

        filename = os.path.join(self.path, self.filename)
        data = pyfits.getdata(filename)
        header = pyfits.getheader(filename)

        w = np.arange(header['naxis3'])
        w = header['CRVAL3'] + header['CDELT3'] * (w - header['CRPIX3'])

        # Convert wavelength from nm to Angstrom
        if w.all() < 1000.:
            w = w * 10

        W = np.reshape(w, (w.size,1,1))
        W = np.repeat(W, data.shape[1], axis=1)
        W = np.repeat(W, data.shape[2], axis=2)

        mean = data.mean(axis=0)
        std = data.std(axis=0)
        mean_mean = np.mean(mean)
        mean_std = np.mean(std)

        # Creating empty space
        mean = - np.inf * np.ones_like(mean)
        std = - np.inf * np.ones_like(mean)
        center = - np.inf * np.ones_like(mean)

        # Filling empty space
        mean = np.where(mean > (mean_mean + mean_std), data.mean(axis=0), -100000)
        std = np.where(mean > (mean_mean + mean_std), data.std(axis=0), -100000)
        center = np.where(mean > (mean_mean + mean_std), ((data * W).sum(axis=0) / data.sum(axis=0)), -100000)

        # fwhm = np.sqrt((data * W * W).sum(axis=0) / data.sum(axis=0) - center ** 2) / (2 * np.sqrt(2 * np.log(2)))

        # p = [0, 0, 0]
        # velocity_map = np.zeros_like(data[0])
        # for i in range(data.shape[2]):
        #     for j in range(data.shape[1]):
        #         print(i, j)
        #
        #         p[0] = data[:,j,i].max()
        #         p[1] = w[np.argmax(data[:,j,i])]
        #         p[2] = fwhm[j,i]
        #
        #         solp, ier = leastsq(self.error_func, p, args=(w,data[:,j,i]))
        #
        #         if ier < 5:
        #             temp = p[1]
        #             temp = (temp - 6562.8) / 6562.8 * 299279.
        #             velocity_map[j,i] = temp
        #         else:
        #             velocity_map[j,i] = -np.inf
        #
        # velocity_map = velocity_map - np.median(velocity_map)

        center = (center - 6562.8) / 6562.8 * 299279.
        # fwhm = (fwhm - 6562.8) / 6562.8 * 299279.

        del header['CRVAL3']
        del header['CDELT3']
        del header['CRPIX3']
        del header['CTYPE3']
        del header['CUNIT3']
        del header['C3_3']
        del header['CD3_3']

        path, filename = os.path.split(filename)
        # pyfits.writeto(os.path.join(path, "vmap_" + filename), velocity_map, header)
        pyfits.writeto(self.safesave(os.path.join(path, "center_" + filename)), center, header)
        pyfits.writeto(self.safesave(os.path.join(path, "std_" + filename)), std, header)
        pyfits.writeto(self.safesave(os.path.join(path, "mean_" + filename)), mean, header)

        return

    def lorentzian(self, p, x):
        return p[0] * (1 / (1 + ((x - p[1]) / p[2])**2))


    def error_func(self, p, x, y):
        return y - self.lorentzian(p, x)

    def safesave(self, name, overwrite=None, verbose=False):
        """
        This is a generic method used to check if a file called 'name' already
        exists. If so, it starts some interaction with the user.

        @param name: the name of the file that will be written in the future.

        @keyword overwrite: if False, this method will interact with the user to
        ask if 'name' file shall be overwritten or if a new name will be given.
        If True, 'name' file is automatically overwritten.

        @keyword verbose: force verbose mode on even when overwrite is automatic.

        v1.0.1 - added 'overwrite' keyword.
               - added 'verbose' keyword.
        """
        import os
        import sys

        v = False if (overwrite is True) else True
        if v:
            print("\n Writing to output file %s\n" % name)

        while os.path.exists(name):

            if overwrite in ['y', 'Y', True]:
                if v or verbose:
                    print(" Overwriting %s file." % name)
                os.remove(name)

            elif overwrite in ['', 'n', 'N', False]:
                name = raw_input("   Please, enter a new filename:\n   > ")

            elif overwrite in ['q']:
                if v:
                    print(" Exiting program.")
                sys.exit()

            else:
                overwrite = raw_input("\tFile\n\t'%s'\n\texist."% name +
                                                    " Overwrite? (y/[n])")
                if v:
                    print(" Writing data-cube to %s" % name)

        return name


if __name__ == '__main__':

    path = '.'
    # filename = 'cube30Dor_SOAR.fits'
    filename = '/home/bquint/public_html/Data/NGC2440.fits'

    main = Main(path, filename)
    main.run()




