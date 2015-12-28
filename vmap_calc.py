#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import print_function, division
import astropy.io.fits as pyfits
import itertools
import logging
import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import leastsq
from scipy import ndimage

__author__ = 'b1quint'
__version__ = '20150406a'

log_formatter = logging.Formatter(" %(levelname)s - %(message)s")

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_formatter)

log = logging.getLogger()
log.addHandler(stream_handler)
log.setLevel(logging.DEBUG)

class Main():

    def __init__(self, filename):
        self.filename = filename
        return

    def fit_continuum(self, data, mask):
        """
        Use the data and the mask to estimate and fit the continuum levels.

        :param data: data-cube in 3D array.
        :param mask: the mask containing the regions where the continuum is dominant.
        :return: array containing the fitted polynomium.
        """
        data = data[:,mask]

        plt.plot(data)
        plt.show()


    def fit_gaussian(self, x, y):
        """
        Fit a gaussian profile to the curve.
        :param x: 1D array with X values.
        :param y: 1D array with Y values.
        :return p: gaussian parameters.
        """
        p = [0, 0, 0, 0]
        p[0] = y.max()
        p[1] = x[y.argmax()]
        p[2] = np.sqrt((y * x * x).sum() / y.sum(axis=0) - p[1] ** 2)
        p[3] = np.median(y)

        err_function = lambda p, x, y: y - self.gaussian(p, x)
        p, s = leastsq(err_function, p, args=(x, y))

        return p

    def fit_lorentzian(self, x, y):
        """
        Fit a lorentzian profile to the curve.
        :param x: 1D array with X values.
        :param y: 1D array with Y values.
        :return p: lorentzian parameters.
        """
        p = [0, 0, 0, 0]
        p[0] = y.max()
        p[1] = x[y.argmax()]
        p[2] = np.sqrt((y * x * x).sum() / y.sum(axis=0) - p[1] ** 2)
        p[3] = np.median(y)

        err_function = lambda p, x, y: y - self.lorentzian(p, x)
        p, s = leastsq(err_function, p, args=(x, y))

        return p

    @staticmethod
    def gaussian(p, x):
        """
        Evaluates a gaussian function.

        :param p: gaussian parameters - p[0] is the amplitude, p[1] is the
        gaussian center, p[2] is related to the FWHM and p[3] is the
        background level.
        :param x: values that will be used to evaluate the gaussian.
        :return: an array with the evaluated gaussian.
        """
        a = p[0]
        b = p[1]
        c = p[2]
        d = p[3]
        g = a * np.exp(-(b - x) **2 / (2*c**2)) + d
        return g

    @staticmethod
    def lorentzian(p, x):
        """
        Evaluates a lorentzian function.

        :param p: lorentzian parameters.
        :param x: values that will be used to evaluate the lorentzian.
        :return: an array with the evaluated lorentzian.
        """
        a = p[0]
        b = p[1]
        c = p[2]
        d = p[3]
        l = a * (1 / (1 + ((x - b) / c)**2)) + d
        return l

    def run(self, show=False):

        # Load data
        d = pyfits.getdata(filename)
        h = pyfits.getheader(filename)

        # Get wavelength calibration
        z = np.arange(h['naxis3'])
        w = h['CRVAL3'] + h['CDELT3'] * (z - h['CRPIX3'])

        # Convert wavelength from nm to Angstrom
        if w.all() < 1000.:
            w *= 10

        # Signal-to-noise clipping
        s = d.sum(axis=2)
        s = s.sum(axis=1)

        gauss_p = self.fit_gaussian(z, s)
        log.debug("Gaussian parameters")
        log.debug("p[0] = %.2f" % gauss_p[0])
        log.debug("p[1] = %.2f" % gauss_p[1])
        log.debug("p[2] = %.2f" % gauss_p[2])
        log.debug("p[3] = %.2f" % gauss_p[3])

        lor_p = self.fit_lorentzian(z, s)
        log.debug("Lorentzian parameters")
        log.debug("p[0] = %.2f" % lor_p[0])
        log.debug("p[1] = %.2f" % lor_p[1])
        log.debug("p[2] = %.2f" % lor_p[2])
        log.debug("p[3] = %.2f" % lor_p[3])

        fwhm = np.abs(gauss_p[2] * 2 * np.sqrt(2 * np.log(2)))
        filter_ = np.where(np.abs(z - gauss_p[1]) < fwhm, True, False)

        if show:
            plt.plot(z, self.gaussian(gauss_p, z), 'r-', lw=2, label='Gaussian Fit')
            plt.plot(z, self.lorentzian(lor_p, z), 'b-', lw=2, label='Lorentzian Fit')
            plt.plot(z, s, 'ko')
            plt.plot(z[filter_], s[filter_], 'ro')
            plt.title('Cube collapsed in XY and fits.')
            plt.grid()
            plt.legend(loc='best')
            plt.show()

        signal = d[filter_].mean(axis=0)
        noise = d[np.logical_not(filter_)].mean(axis=0)
        snr = signal / noise
        snr = ndimage.median_filter(snr, 3)
        snr_laplacian = ndimage.morphological_laplace(snr, size=3)

        snr_mask = np.where(snr > 3, True, False)
        snr_mask *= np.where(np.abs(snr_laplacian) < 0.90, True, False)
        snr_mask = ndimage.binary_erosion(snr_mask)
        snr_mask = ndimage.binary_erosion(snr_mask)
        snr_mask = ndimage.binary_dilation(snr_mask)
        snr_mask = ndimage.binary_dilation(snr_mask)

        pyfits.writeto(filename.replace('.','.SNR.'), snr, h, clobber=True)
        pyfits.writeto(filename.replace('.','.SNR_LAPLACIAN.'), snr_laplacian, h, clobber=True)

        if show:

            fig1 = plt.figure(figsize=(15,5))

            plt.title('Signal-to-Noise Ratio')

            ax1 = plt.subplot2grid((1, 3),(0, 0))
            im1 = ax1.imshow(snr, cmap='cubehelix_r', interpolation='nearest',
                       origin='lower', vmin=2.5, vmax=25)
            ax1.set_title('Raw')
            cbar1 = plt.colorbar(mappable=im1, ax=ax1, use_gridspec=True, orientation='vertical')

            ax2 = plt.subplot2grid((1, 3),(0, 1))
            im2 = ax2.imshow(snr_mask, cmap='gray_r', interpolation='nearest',
                       origin='lower', vmin=0, vmax=1)
            ax2.set_title('Mask')
            cbar2 = plt.colorbar(mappable=im2, ax=ax2, use_gridspec=True, orientation='vertical')

            ax3 = plt.subplot2grid((1, 3),(0, 2))
            im3 = ax3.imshow(snr_mask * snr, cmap='cubehelix_r', interpolation='nearest',
                       origin='lower', vmin=0, vmax=25)
            ax3.set_title('Masked')
            cbar3 = plt.colorbar(mappable=im3, ax=ax3, use_gridspec=True, orientation='vertical')

            plt.tight_layout()
            plt.show()

        # Adjust continuum
        continuum = self.fit_continuum(d, snr_mask)

        # flux = np.zeros_like(snr)
        # velocity = np.zeros_like(snr)
        # width = np.zeros_like(snr)

        # for (j, i) in itertools.product(range(d.shape[1]), range(d.shape[2])):
        #
        #     if not snr_mask[j, i]:







        # W = np.reshape(w, (w.size,1,1))
        # W = np.repeat(W, data.shape[1], axis=1)
        # W = np.repeat(W, data.shape[2], axis=2)
        #
        # mean = data.mean(axis=0)
        # std = data.std(axis=0)
        # mean_mean = np.mean(mean)
        # mean_std = np.mean(std)
        #
        # # Creating empty space
        # mean = - np.inf * np.ones_like(mean)
        # std = - np.inf * np.ones_like(mean)
        # center = - np.inf * np.ones_like(mean)
        #
        # # Filling empty space
        # mean = data.mean(axis=0)
        # std = data.std(axis=0)
        # center = ((data * W).sum(axis=0) / data.sum(axis=0))
        # fwhm = np.sqrt((data * W * W).sum(axis=0) / data.sum(axis=0) - center ** 2) / (2 * np.sqrt(2 * np.log(2)))
        #
        # # Getting SNR
        # signal = data[10:25].sum(axis=0)
        # noise = data[25:].sum(axis=0)
        # snr = signal / noise
        #
        # # Correcting images
        # target_snr = 10
        # mean = np.where(snr > target_snr, mean, -10000)
        # std = np.where(snr > target_snr, std, -10000)
        # center = np.where(snr > target_snr, center, -10000)
        # fwhm = np.where(snr > target_snr, fwhm, -10000)


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

        # center = (center - 6562.8) / 6562.8 * 299279.
        # fwhm = (fwhm - 6562.8) / 6562.8 * 299279.

        # del header['CRVAL3']
        # del header['CDELT3']
        # del header['CRPIX3']
        # del header['CTYPE3']
        # del header['CUNIT3']
        # del header['C3_3']
        # del header['CD3_3']
        #
        # path, filename = os.path.split(filename)
        # # pyfits.writeto(os.path.join(path, "vmap_" + filename), velocity_map, header)
        # pyfits.writeto(self.safesave(os.path.join(path, "center_" + filename)), center, header, clobber=True)
        # pyfits.writeto(self.safesave(os.path.join(path, "std_" + filename)), std, header, clobber=True)
        # pyfits.writeto(self.safesave(os.path.join(path, "mean_" + filename)), mean, header, clobber=True)
        # pyfits.writeto(self.safesave(os.path.join(path, "fwhm_" + filename)), fwhm, header, clobber=True)
        #
        return





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

    filename = '/home/b1quint/Data/NGC2440.fits'
    main = Main(filename)
    main.run(show=True)




