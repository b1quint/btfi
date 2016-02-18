#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import print_function, division

from scipy.signal.cont2discrete import cont2discrete

import astropy.io.fits as pyfits
import itertools
import logging
import matplotlib.pyplot as plt
import numpy as np
import os

from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy.optimize import leastsq
from scipy import ndimage

from matplotlib import colors
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Rectangle
from matplotlib import widgets
from mpl_toolkits.axes_grid1 import make_axes_locatable

__author__ = 'b1quint'
__version__ = '20150406a'

class Main():

    def __init__(self, filename):

        self.filename = filename
        self.data = None
        self.header = None
        self.rect = None

        self.x1 = 0
        self.x2 = 0
        self.y1 = 0
        self.y2 = 0

        self.ax1 = None
        self.ax2 = None
        self.RS = None

        return

    @staticmethod
    def double_gaussian(p, x):
        """
        Evaluates a double gaussian function.

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
        e = p[4]
        f = p[5]
        g = p[6]
        gauss = a * np.exp(-(b - x) **2 / (2*c**2)) + \
                d * np.exp(-(e - x) **2 / (2*f**2)) + g
        return gauss

    def fit_continuum(self, data):
        """
        Use the data and the mask to estimate and fit the continuum levels.

        :param collapsed_data: data-cube in 3D array.
        :return: array containing the fitted polynomium.
        """
        collapsed_data = np.mean(data, axis=0)
        norm = ImageNormalize(vmin=collapsed_data.min(), vmax=collapsed_data.max(), stretch=LogStretch())

        fig = plt.figure(figsize=(8, 6))
        fig.suptitle('Draw a rectangle using the mouse. \nPress <ENTER> to ' +
                      'accept it and carry on or "Q" to leave.')

        gs = GridSpec(3, 3, height_ratios=[1, 12, 3], width_ratios=[1, 10, 1])
        ax1 = plt.subplot(gs[4])
        im1 = ax1.imshow(collapsed_data, origin='lower', interpolation='nearest',
                         cmap='hot_r', norm=norm)
        ax1.grid()

        self.ax2 = plt.subplot(gs[7])
        self.ax2.xaxis.set_ticklabels([])
        self.ax2.yaxis.set_ticklabels([])

        self.RS = MyRectangleSelector(ax1, self.line_select_callback, drawtype='box',
                                 useblit=True, button=[1], minspanx=5,
                                 minspany=5, spancoords='pixels',
                                 rectprops = dict(facecolor='green',
                                                  edgecolor = 'green',
                                                  alpha=0.5, fill=True))


        self.RS.set_active(True)
        self.RS.connect_event('button_press_event', self.on_mouse_click)
        self.RS.connect_event('button_release_event', lambda e: self.on_mouse_release(e, data))
        self.RS.connect_event('key_press_event', self.on_key_press)

        gs.tight_layout(fig)
        plt.show()

        x1 = min(self.x1, self.x2)
        x2 = max(self.x1, self.x2)
        y1 = min(self.y1, self.y2)
        y2 = max(self.y1, self.y2)

        if x1 == x2:
            log.warning('x1 and x2 are the same. Using the whole image width.')
            x1 = 0
            x2 = -1

        if y1 == y2:
            log.warning('y1 and y2 are the same. Using the whole image height.')
            y1 = 0
            y2 = -1

        data = data[:, y1:y2, x1:x2]
        data = data.mean(axis=1)
        data = data.mean(axis=1)
        median, std = np.median(data), np.std(data)
        c = np.where(np.abs(data - median) < std, True, False)
        x = np.arange(data.size)
        p = np.polyfit(x[c], data[c], 3)
        y = np.polyval(p, x)

        return y


    def fit_double_gaussian(self, x, y):
        """
        Fit a double gaussian profile to the curve.
        :param x: 1D array with X values.
        :param y: 1D array with Y values.
        :return p: gaussian parameters.
        """
        p = [1, 1, 1, 1, 1, 1, 1]

        def err_function(p, x, y):
            return y - self.double_gaussian(p, x)

        p, s = leastsq(err_function, p, args=(x, y))
        return p, s


    def fit_gaussian(self, x, y):
        """
        Fit a gaussian profile to the curve.
        :param x: 1D array with X values.
        :param y: 1D array with Y values.
        :return p: gaussian parameters.
        """
        p = [0, 0, 0, 0]
        p[0] = y.max()
        p[1] = x[x.size // 2]
        p[2] = 1
        p[3] = np.median(y)

        def err_function(p, x, y):
            return y - self.gaussian(p, x)

        p, s = leastsq(err_function, p, args=(x, y))
        return p, s

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

        def err_function(p, x, y):
            return y - self.lorentzian(p, x)

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

    def line_select_callback(self, eclick, erelease):
        'eclick and erelease are the press and release events'
        log.debug("Click at [{:},{:}]. Release at [{:},{:}]".format(
            eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata
        ))
        self.x1 = int(eclick.xdata)
        self.y1 = int(eclick.ydata)
        self.x2 = int(erelease.xdata)
        self.y2 = int(erelease.ydata)
        return

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

    def on_key_press(self, event):

        log.debug('Key pressed: ' + event.key)
        if event.key.lower() in ['enter', ' ']:
            curr_fig = plt.gcf()
            plt.close(curr_fig)
        if event.key.lower() in ['q', 'escape']:
            from sys import exit
            exit()

        return

    def on_mouse_click(self, event):

        try:
            self.rect.remove()
        except (NotImplementedError, AttributeError):
            log.debug('No existing rectangle. Creating one now.')

        curr_fig = plt.gcf()
        curr_fig.canvas.draw()

        self.x1 = int(event.xdata)
        self.y1 = int(event.ydata)

        log.debug('Mouse press event at [%d, %d]' % (self.x1, self.y1))
        log.debug(event.inaxes)
        log.debug(event.xdata)
        log.debug(event.ydata)

        return

    def on_mouse_release(self, event, data):
        """
        Event handler to work while selecting a continuum area.

        @param event: 'mouse_release' event
        @param data: the 3D datacube to be collapsed in XY directions.
        @return: None
        """
        x = int(event.xdata)
        y = int(event.ydata)

        log.debug('Mouse release event at [%d, %d]' % (x, y))

        self.rect = Rectangle((self.x1, self.y1),
                              self.x2 - self.x1,
                              self.y2 - self.y1,
                              facecolor='blue',
                              edgecolor='blue',
                              alpha=0.3)

        log.debug(event.inaxes)
        log.debug(event.xdata)
        log.debug(event.ydata)
        log.debug(self.x1)
        log.debug(self.x2)
        log.debug(self.y1)
        log.debug(self.y2)

        x1 = min(self.x1, self.x2)
        x2 = max(self.x1, self.x2)
        y1 = min(self.y1, self.y1)
        y2 = max(self.y1, self.y2)

        data = data[:, y1:y2, x1:x2]
        data = data.mean(axis=1)
        data = data.mean(axis=1)
        median, std = np.median(data), np.std(data)
        c = np.where(np.abs(data - median) < std, True, False)
        x = np.arange(data.size)
        p = np.polyfit(x[c], data[c], 3)
        fx = np.linspace(0, data.size, data.size * 4)
        fy = np.polyval(p, fx)

        self.ax2.cla()
        self.ax2.plot(x[c], data[c], 'ko')
        self.ax2.plot(x[~c], data[~c], 'ro', alpha=0.25)
        self.ax2.plot(fx, fy, 'r-')
        self.ax2.plot(fx, np.ones_like(fx), 'r--')
        self.ax2.plot(fx, np.ones_like(fx) * median - std, 'r--', alpha=0.25)
        self.ax2.plot(fx, np.ones_like(fx) * median + std, 'r--', alpha=0.25)
        self.ax2.set_xlim(x.min(), x.max())
        self.ax2.locator_params(axis='y', nbins=2)

        curr_ax = event.inaxes
        curr_fig = plt.gcf()
        curr_ax.add_patch(self.rect)
        curr_fig.canvas.draw()

        return

    def run(self, show=False):

        # Load data
        filename = self.filename
        d = pyfits.getdata(filename)
        h = pyfits.getheader(filename)
        path, filename = os.path.split(filename)

        # Get wavelength calibration
        z = np.arange(h['naxis3'])
        w = h['CRVAL3'] + h['CDELT3'] * (z - h['CRPIX3'])

        # Convert wavelength from nm to Angstrom
        # if w.all() < 1000.:
        #     w *= 10

        # Signal-to-noise clipping
        s = d.sum(axis=2)
        s = s.sum(axis=1)

        gauss_p, _ = self.fit_gaussian(z, s)
        log.debug("Gaussian parameters ---")
        log.debug("p[0] = %.2f" % gauss_p[0])
        log.debug("p[1] = %.2f" % gauss_p[1])
        log.debug("p[2] = %.2f" % gauss_p[2])
        log.debug("p[3] = %.2f" % gauss_p[3])

        lor_p = self.fit_lorentzian(z, s)
        log.debug("Lorentzian parameters ---")
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
            plt.gcf().canvas.mpl_connect('key_press_event', self.on_key_press)
            plt.show()

        signal = d[filter_].mean(axis=0)
        noise = d[np.logical_not(filter_)].mean(axis=0)
        snr = signal / noise
        snr = ndimage.median_filter(snr, 3)
        snr_mask = np.where(snr > 2, True, False)

        snr_laplacian = ndimage.morphological_laplace(snr * snr_mask, size=3)
        snr_mask *= np.where(np.abs(snr_laplacian) < 5.0, True, False)

        snr_mask = ndimage.binary_opening(snr_mask, iterations=5)
        snr_mask = ndimage.binary_closing(snr_mask, iterations=5)
        if show:

            fig1 = plt.figure(figsize=(20, 5))
            plt.title('Signal-to-Noise Ratio')
            gs = GridSpec(1, 3)

            ax1 = plt.subplot(gs[0])
            ax1.set_title('SNR')
            im1 = ax1.imshow(snr, cmap='cubehelix', interpolation='nearest',
                       origin='lower', vmin=0)
            div1 = make_axes_locatable(ax1)
            cax1 = div1.append_axes("right", size="5%", pad=0.05)
            cbar1 = plt.colorbar(mappable=im1, cax=cax1, use_gridspec=True, 
                                 orientation='vertical')

            ax2 = plt.subplot(gs[1])
            ax2.set_title('Mask')
            im2 = ax2.imshow(np.where(snr_mask, 1, 0), cmap='gray',
                             interpolation='nearest', origin='lower',
                             vmin=0, vmax=1)
            div2 = make_axes_locatable(ax2)
            cax2 = div2.append_axes("right", size="5%", pad=0.05)
            cbar2 = plt.colorbar(mappable=im2, cax=cax2, use_gridspec=True, 
                                 orientation='vertical')

            cmap = plt.get_cmap('cubehelix')
            cmap.set_bad('w', 1.0)
            ax3 = plt.subplot(gs[2])
            ax3.set_title('Masked')
            im3 = ax3.imshow(np.ma.masked_where(~snr_mask, snr), cmap=cmap, interpolation='nearest',
                       origin='lower', vmin=0)
            div3 = make_axes_locatable(ax3)
            cax3 = div3.append_axes("right", size="5%", pad=0.05)
            cbar3 = plt.colorbar(mappable=im3, cax=cax3, use_gridspec=True, orientation='vertical')

            plt.gcf().canvas.mpl_connect('key_press_event', self.on_key_press)
            gs.tight_layout(fig1)
            plt.show()

        pyfits.writeto(filename.replace('.','.SNR.'), snr, h, clobber=True)
        pyfits.writeto(filename.replace('.','.SNR_LAPLACIAN.'), snr_laplacian, h, clobber=True)

        # Adjust continuum
        continuum = self.fit_continuum(d)

        # Subtract continuum
        continuum = np.reshape(continuum, (continuum.size, 1, 1))
        continuum = np.repeat(continuum, d.shape[1], axis=1)
        continuum = np.repeat(continuum, d.shape[2], axis=2)
        d -= continuum
        del continuum

        # Save the data-cube with the continuum subtracted
        pyfits.writeto(os.path.join(path,
                                    filename.replace('.fits', '.no_cont.fits')),
                       d, h, clobber=True)

        # Adjust each pixel to a gaussian
        flux = np.zeros_like(snr) - np.inf
        velocity = np.zeros_like(snr) - np.inf
        width = np.zeros_like(snr) - np.inf
        continuum = np.zeros_like(snr) - np.inf
        niter = np.zeros_like(snr) - np.inf
        fit_gauss_cube = np.zeros_like(d)

        for (j, i) in itertools.product(range(d.shape[1]), range(d.shape[2])):

            if snr_mask[j, i]:
                p, s = self.fit_gaussian(w, d[:,j,i])
            else:
                p, s = [-np.inf, -np.inf, -np.inf, -np.inf], -1

            flux[j, i] = p[0]
            velocity[j, i] = p[1]
            width[j, i] = p[2]
            continuum[j, i] = p[3]
            niter[j,i] = s

            fit_gauss_cube[:, j, i] = self.gaussian(p, w)

            log.debug("%04d %04d" % (j, i))

        mask = np.ma.masked_equal(snr_mask, False)

        pyfits.writeto(os.path.join(path,
                                    filename.replace('.fits', '.GAUSS.fits')),
                       fit_gauss_cube, h, clobber=True)

        if show:

            fig1 = plt.figure(figsize=(16,6))
            plt.title('Final Results')
            gs = GridSpec(2, 3)

            magma = plt.get_cmap('magma')
            magma.set_bad('w', 1.0)
            ax1 = plt.subplot(gs[0])
            ax1.set_title('Collapsed cube')
            im1 = ax1.imshow(np.ma.masked_where(~snr_mask, d.sum(axis=0)),
                             cmap=magma, interpolation='nearest',
                             origin='lower')

            divider1 = make_axes_locatable(ax1)
            cax1 = divider1.append_axes("right", size="5%", pad=0.05)
            cbar1 = plt.colorbar(mappable=im1, cax=cax1, use_gridspec=True,
                                 orientation='vertical')

            viridis = plt.get_cmap('viridis')
            viridis.set_bad('w', 1.0)
            ax2 = plt.subplot(gs[1])
            ax2.set_title('Fit n-iteractions')
            im2 = ax2.imshow(np.ma.masked_where(~snr_mask, niter),
                             cmap=viridis, interpolation='nearest',
                             origin='lower')

            divider2 = make_axes_locatable(ax2)
            cax2 = divider2.append_axes("right", size="5%", pad=0.05)
            cbar2 = plt.colorbar(mappable=im2, cax=cax2, use_gridspec=True,
                                 orientation='vertical')

            ax3 = plt.subplot(gs[3])
            ax3.set_title('Flux')

            divider3 = make_axes_locatable(ax3)
            cax3 = divider3.append_axes("right", size="5%", pad=0.05)
            im3 = ax3.imshow(flux * mask, cmap=viridis,
                             interpolation='nearest', origin='lower')

            cbar3 = plt.colorbar(mappable=im3, cax=cax3, use_gridspec=True,
                                  orientation='vertical')

            RdYlBu = plt.get_cmap('RdYlBu')
            RdYlBu.set_bad('w', 1.0)

            ax4 = plt.subplot(gs[4])
            ax4.set_title('Center')

            vmin = np.median(velocity * mask) - 2 * np.std(velocity * mask)
            vmax = np.median(velocity * mask) + 2 * np.std(velocity * mask)
            im4 = ax4.imshow(velocity * mask, cmap=RdYlBu,
                             interpolation='nearest', origin='lower',
                             vmin=vmin, vmax=vmax)

            divider4 = make_axes_locatable(ax4)
            cax4 = divider4.append_axes("right", size="5%", pad=0.05)
            cbar4 = plt.colorbar(mappable=im4, cax=cax4, use_gridspec=True,
                                 orientation='vertical')

            ax5 = plt.subplot(gs[5])
            ax5.set_title('Width')

            vmin = np.median(width * mask) - 2 * np.std(width * mask)
            vmax = np.median(width * mask) + 2 * np.std(width * mask)
            im5 = ax5.imshow(width * mask, cmap=viridis,
                             interpolation='nearest', origin='lower')

            divider5 = make_axes_locatable(ax5)
            cax5 = divider5.append_axes("right", size="5%", pad=0.05)
            cbar5 = plt.colorbar(mappable=im5, cax=cax5, use_gridspec=True,
                                 orientation='vertical')


            plt.gcf().canvas.mpl_connect('key_press_event', self.on_key_press)
            plt.tight_layout()
            plt.show()


        del h['CRVAL3']
        del h['CDELT3']
        del h['CRPIX3']
        del h['CTYPE3']
        del h['CUNIT3']
        del h['C3_3']
        del h['CD3_3']
        
        pyfits.writeto(
            os.path.join(path, filename.replace('.fits', '.flux.fits')),
            flux, h, clobber=True)

        pyfits.writeto(
            os.path.join(path, filename.replace('.fits', '.velo.fits')),
            velocity, h, clobber=True)

        pyfits.writeto(
            os.path.join(path, filename.replace('.fits', '.width.fits')),
            width, h, clobber=True)

        pyfits.writeto(
            os.path.join(path, filename.replace('.fits', '.cont.fits')),
            continuum, h, clobber=True)

        return

    @staticmethod
    def safesave(name, overwrite=None, verbose=False):
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

class MyColors:

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def __init__(self):
        pass

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        return

    def enable(self):
        self.HEADER = '\033[95m'
        self.OKBLUE = '\033[94m'
        self.OKGREEN = '\033[92m'
        self.WARNING = '\033[93m'
        self.FAIL = '\033[91m'
        self.ENDC = '\033[0m'
        return


class MyFormatter(logging.Formatter):
    """
    Custom log format for a better interaction with the user.
    """
    info_fmt = u"    %(msg)s"

    wrn_fmt = MyColors.WARNING + "[W]" + MyColors.ENDC + \
        u" %(msg)s"

    dbg_fmt = u"{0}[D]{1} %(asctime)s %(filename)s %(funcName)s"  + \
              u"  %(lineno)d: %(msg)s".format(
        MyColors.OKBLUE, MyColors.ENDC)

    err_fmt = u"{0}[E]{1} %(asctime)s %(filename)s %(funcName)s" + \
              u"  %(lineno)d: %(msg)s".format(
        MyColors.FAIL, MyColors.ENDC)

    def __init__(self, fmt=u"%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)
        return


    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = MyFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = MyFormatter.err_fmt

        elif record.levelno == logging.WARNING:
            self._fmt = MyFormatter.wrn_fmt

        elif record.levelno == logging.CRITICAL:
            self._fmt = MyFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


class MyRectangleSelector(widgets.RectangleSelector):
    """
    Custom rectangle selector used to make keep the selected area
    persistent on the screen. Based on:

    http://stackoverflow.com/a/34517699/2333908
    """
    def release(self, event):
        super(MyRectangleSelector, self).release(event)
        self.to_draw.set_visible(True)
        self.canvas.draw()

# Main thread ---
log_formatter = logging.Formatter(" %(levelname)s - %(message)s")
log_formatter = MyFormatter()

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_formatter)

log = logging.getLogger()
log.addHandler(stream_handler)
log.setLevel(logging.DEBUG)

if __name__ == '__main__':

    # filename = '/home/b1quint/Data/NGC2440.fits' # lemonjuice
    filename = '/home/bquint/Data/NGC2440_x200x700_y400y700.fits' # soarbr3
    main = Main(filename)
    main.run(show=True)




