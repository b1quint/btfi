#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    Fabry-Perot Wavelength Calibration

    This file needs at least two identified emission lines extracted from one
    or more fits files containing data-cubes or spectra.

    by Bruno C. Quint <bquint at ctio dot noao dot edu>
    2015.09.11
"""
from __future__ import division
import argparse
import astropy.io.fits as pyfits
import logging
import matplotlib.pyplot as plt
import numpy as np
import sys

from scipy.signal import argrelextrema
from scipy.interpolate import splev, splrep

def main():

    # Create and setup log ---
    fmt = MyFormatter()
    hdlr = logging.StreamHandler(sys.stdout)
    hdlr.setFormatter(fmt)
    logging.root.addHandler(hdlr)
    logging.root.setLevel(logging.DEBUG)
    log = logging.getLogger(__name__)

    # Parse command line arguments ---
    parser = argparse.ArgumentParser(
        description="Calculate the wavelength calibration for FP observation.")
    parser.add_argument('-v', '--verbose_level', default=1, type=int,
                        help="Set verbose level: 0 - Quiet, 1 - Info," + \
                             "2 - Warnings, 3 - Debug. (Default: 1)")
    parser.add_argument('input_files', type=str, nargs='+',
                        help="Input data-cube(s) for calibration.")

    args = parser.parse_args()

    # Defining verbose/logging level ---
    if args.verbose_level == 0:
        log.setLevel(logging.CRITICAL)
    elif args.verbose_level == 1:
        log.setLevel(logging.INFO)
    elif args.verbose_level == 2:
        log.setLevel(logging.WARNING)
    elif args.verbose_level == 3:
        log.setLevel(logging.DEBUG)
    else:
        log.error("Invalid log level. Setting DEBUG level.")
        log.setLevel(logging.DEBUG)

    log.info("\n    Fabry-Perot Wavelength Calibration\n" + \
             "    by Bruno C. Quint - bquint at ctio dot noao dot edu\n" + \
             "    Version 0.0a\n")

    # Read files and extract spectra ---
    if len(args.input_files) == 0:
        log.error("You must set at leat one file. Leaving now.")
        sys.exit()

    list_of_spectra = []
    for input_file in args.input_files:
        log.info("Reading file {}".format(input_file))
        try:
            h = pyfits.getheader(input_file)
        except IOError:
            log.error("{} not found.".format(input_file))
            continue

        if h['NAXIS'] == 1:
            s = pyfits.getdata(input_file)
            try:
                z = (np.arange(h['NAXIS1']) - h['CRPIX1'] + 1) * \
                    h['CDELT1'] + h['CRVAL1']
            except KeyError:
                log.warning("Spectral calibration not found in " + \
                            "{}'s header.".format(input_file))
                z = np.arange(h['NAXIS1'])
        else:
            log.warning()

        # Smooth spectrum --
        s = np.convolve(s, [0.1, 0.15, 0.5, 0.15, 0.1], mode='same')

        # Interpolate --
        tck = splrep(z, s, k=5)
        znew = np.linspace(z[0], z[-1], num=1000, endpoint=True)
        snew = splev(znew, tck, der=0)

        # Plot --
        log.info("A new window has opened. Please, click near the peak(s)" + \
                 " that you want to identify.")

        # plt.ion()
        p = PointPicker(z, s, znew, snew)
        plt.show()

    return 0


class PointPicker(object):

    log = logging.getLogger(__name__)

    def __init__(self,x, y, xnew, ynew):

        self.log.debug("Creating object.")
        self.x = xnew
        self.y = ynew

        ptp = np.ptp(ynew)
        ymin = np.min(ynew) - 0.10 * ptp
        ymax = np.max(ynew) + 0.10 * ptp

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        self.points2d, = self.ax.plot(x, y, 'bo')
        self.lines2d, = self.ax.plot(xnew, ynew, 'b-', picker=10)

        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('key_press_event', self.onpress)

        self.ax.set_ylim(ymin=ymin, ymax=ymax)
        self.ax.grid()
        self.log.debug("Fininsh creation and leaving constructor.")
        return


    def onpress(self, event):
        """define some key press events"""
        if event.key.lower() == 'q':
            sys.exit()

        if event.key.lower() == 'c':
             self.clear_lines()

        if event.key.lower() == 'd':
             self.delete_line()

        if event.key.lower() == '':
            self.accept_lines()

        if event.key.lower() == '':
            self.refuse_spectrum()


    def onpick(self,event):
        """Define what happens when one clicks near the data"""
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        self.log.debug('Object picked at x = {}, y = {}'.format(x, y))

        i = np.argmin(np.abs(x - self.x))
        self.log.debug('Nearest index: {}'.format(i))

        i0 = max(0, i-2)
        iz = min(i+3, len(self.x))

        slice_index = np.arange(0, 5) - 2
        index = -2
        while index != 0:

            slice_index = slice_index[i + slice_index >= 0] # Protect right border
            slice_index = slice_index[i + slice_index <= len(self.x)] # Protect left border
            print(i + slice_index)
            yy = self.y[i + slice_index]
            index = np.argmax(yy) - 2
            self.log.debug('Finding max index: {}, {}'.format(i, index))

            i = i + 1 if index > 0 else i - 1

        L =  self.ax.axvline(x=self.x[i], c='red')
        self.fig.canvas.draw()


class TColors:

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

# Custom formatter
class MyFormatter(logging.Formatter):

    info_fmt = "    %(msg)s"
    wrn_fmt  = TColors.WARNING + "[W]" + TColors.ENDC + \
               " %(msg)s"
    dbg_fmt  = TColors.OKBLUE + "[D]" + TColors.ENDC + \
               " %(asctime)s %(filename)s %(funcName)s " + \
               " %(lineno)d: %(msg)s"
    err_fmt  = TColors.FAIL + "[E]" + TColors.ENDC + \
               " %(asctime)s %(filename)s %(funcName)s " + \
               " %(lineno)d: %(msg)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
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
    main()