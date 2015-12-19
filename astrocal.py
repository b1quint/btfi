#/usr/bin/python
# -*- coding: utf8 -*-
"""
    Interactive Astrometric Calibration
"""
from __future__ import division, print_function

import aplpy
import argparse
import astropy.io.fits as pyfits
import logging
import matplotlib.pyplot as plt
import numpy as np

from astropy import coordinates
from astropy import units as u
from matplotlib.patches import Rectangle
from subprocess import call

__author__ = 'Bruno Quint'

class Main(object):

    plate_scale = 0.0455 * u.arcsec
    _2MASS = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'+ \
               'outfmt=1&objstr={:5f}+{:5f}&spatial=Cone&radius=400&' + \
               'catalog=fp_psc'

    def __init__(self):
        return

    def run(self, args):

        # Open File
        filename = args.filename
        data = pyfits.getdata(filename)
        header = pyfits.getheader(filename)

        # Show file
        self.ax1.imshow(data, interpolation='nearest', origin='lower', cmap='gray_r')
        self.ax1.xaxis.set_ticklabels([])
        self.ax1.yaxis.set_ticklabels([])

        # Get object coordinated from header
        ra, dec = header['OBSRA'], header['OBSDEC']
        c = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))

        # Download cathalog
        # call(['wget',
        #       self._2MASS.format(c.ra.deg, c.dec.deg),
        #       '-O',
        #       '_2MASS.tbl'])
        ra, dec, j_m = np.loadtxt('_2MASS.tbl', skiprows=103, usecols=[0, 1, 8], unpack=True)

        # Showing detector area
        width = np.ones(ra.size) * 3. / 60. * u.deg
        height = np.ones(ra.size) * 3. / 60. * u.deg
        self.ax2.scatter(ra, dec, s=j_m, alpha=0.5, facecolor='None', edgecolor='red')
        self.ax2.grid()

        plt.show()

        return


# Terminal Colors ---
class TColors:

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


# Custom formatter ---
class MyFormatter(logging.Formatter):

    info_fmt = u"    %(msg)s"
    wrn_fmt = TColors.WARNING + "[W]" + TColors.ENDC + \
        u" %(msg)s"
    dbg_fmt = u"{0}[D]{1} %(asctime)s %(filename)s %(funcName)s" + \
        u"%(lineno)d: %(msg)s".format(TColors.OKBLUE, TColors.ENDC)
    err_fmt = u"{0}[E]{1} %(asctime)s %(filename)s %(funcName)s  %(lineno)d:" +\
        u"%(msg)s".format(TColors.FAIL, TColors.ENDC)

    def __init__(self, fmt=u"%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

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


if __name__ == '__main__':

   # Parsing Arguments ---
    parser = argparse.ArgumentParser(
        description="Interactive way of finding an astrometric calibration " +
                    "solution for FITS images.")

    parser.add_argument('filename', type=str,
                        help="input filename.")

    parser.add_argument('-b', '--binning', type=int, default=1,
                        help='Image binning.')

    parser.add_argument('-v','--verbose', type=int, default=0,
                        help="Verbose level: 0 - None, 1 - Info, 2 - Debug.")

    args = parser.parse_args()

    # Main Thread ---
    main = Main()
    main.run(args)
