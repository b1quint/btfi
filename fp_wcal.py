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
# from builtins import input

import argparse
import astropy.io.fits as pyfits
import logging
import matplotlib.pyplot as plt
import numpy as np
import sys

from scipy.interpolate import splev, splrep

import matplotlib as mpl
mpl.rcParams['toolbar'] = 'None'


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
        description=u'Calculate the wavelength calibration for FP observation.')
    parser.add_argument('-v', '--verbose_level', default=1, type=int,
                        help=u'Set verbose level: 0 - Quiet, 1 - Info,' +
                             u"2 - Warnings, 3 - Debug. (Default: 1)")
    parser.add_argument("input_files", type=str, nargs='+',
                        help=u"Input data-cube(s) for calibration.")

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
        log.error(u"Invalid log level. Setting DEBUG level.")
        log.setLevel(logging.DEBUG)

    log.info(u"\n    Fabry-Perot Wavelength Calibration\n" + \
             u"    by Bruno C. Quint - bquint at ctio dot noao dot edu\n" + \
             u"    Version 0.0a\n")

    # Read files and extract spectra ---
    if len(args.input_files) == 0:
        log.error(u"You must set at leat one file. Leaving now.")
        sys.exit()

    list_of_spectra = []
    pairs = {}
    for input_file in args.input_files:
        log.info(u"Reading file {}".format(input_file))
        try:
            h = pyfits.getheader(input_file)
        except IOError:
            log.error(u"{} not found.".format(input_file))
            continue

        if h['NAXIS'] == 1:
            s = pyfits.getdata(input_file)
            try:
                z = (np.arange(h['NAXIS1']) - h['CRPIX1'] + 1) * \
                    h['CDELT1'] + h['CRVAL1']
            except KeyError:
                log.warning(u"Spectral calibration not found in " + \
                            u"{}'s header.".format(input_file))
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
        log.info(u"A new window has opened. Please, click near the peak(s)" + \
                 " that you want to identify.")
        log.info(u"Click near a peak to add it.")
        log.info(u"Press:")
        log.info(u"<space_bar> or <ENTER> to accept the lines selected.")
        log.info(u'<d> near a line to delete it.')
        log.info(u'<c> to clear all lines and start over.')
        log.info(u'<ESCAPE> to ignore the current file and jump to the next.')
        log.info(u'<Q> to leave the program.')

        # plt.ion()
        p = PointPicker(input_file, z, s, znew, snew)

        for x in p.xs:
            w = input(u"Type a wavelength for the line at " +
                           u"{:0.1f}:\n> ".format(x))
            if w in pairs.keys():
                pairs[w].append(x)
            else:
                pairs[w] = [x]

    # Calculating the wavelength calibration ---
    log.debug(pairs)
    curvature = 0
    while curvature not in [1,2]:
        curvature = raw_input(u"    Do the rings increase or decrease with Z?" +
            u"    1 - Increase\n    2 - Decrease\n > ")

    return 0


class PointPicker(object):
    log = logging.getLogger(__name__)

    def __init__(self, filename, x, y, xnew, ynew):

        self.log.debug(u"Creating object.")
        self.x = xnew
        self.y = ynew

        ptp = np.ptp(ynew)
        ymin = np.min(ynew) - 0.10 * ptp
        ymax = np.max(ynew) + 0.10 * ptp

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        self.points2d, = self.ax.plot(x, y, 'bo')
        self.lines2d, = self.ax.plot(xnew, ynew, 'b-', picker=10)
        self.vlines = []
        self.xs = []

        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('key_press_event', self.onpress)

        self.ax.set_title(filename)
        self.ax.set_ylim(ymin=ymin, ymax=ymax)
        self.ax.set_xlabel(u"Spectral Units [?]")
        self.ax.set_ylabel(u"Intensity [?]")
        self.ax.grid()
        self.log.debug(u"Fininsh creation and leaving constructor.")

        plt.show()
        return

    def accept_lines(self, event):
        """Accept the current lines drawn"""
        for line in self.ax.lines[2:]:
            x, _ = line.get_xdata()
            self.xs.append(x)
        plt.close('all')
        return

    def clear_lines(self, event):
        """Delete all vertical lines in the plot."""
        for l in self.vlines:
            try:
                l.remove()
            except ValueError as e:
                self.log.error(u"Runtime error - " + e.message)
                pass
        self.fig.canvas.draw()
        return

    def delete_line(self, event):
        """Delete a single line near the cursor."""
        for line in self.ax.lines[2:]:
            xdata, _ = line.get_xdata()
            self.log.debug(u"xdata = {0:0.2f}".format(xdata) +
                           u", event.xdata = {0:0.2f}".format(event.xdata))
            xdata_index = np.argmin(np.abs(self.x - xdata))
            exdata_index = np.argmin(np.abs(self.x - event.xdata))
            self.log.debug(u"xdata_i = {:d}".format(xdata_index) +
                           u", event.xdata_i = {:d}".format(exdata_index))
            if np.abs(xdata_index - exdata_index) < 10:
                try:
                    line.remove()
                except ValueError:
                    pass
            self.fig.canvas.draw()
        return

    def onpress(self, event):
        """define some key press events"""
        self.log.debug(u"Key pressed: " + event.key.lower())
        if event.key == 'Q':
            self.log.info(u'You pressed "Q". Leaving the program now.')
            self.log.info(u'Bye!')
            sys.exit()

        if event.key.lower() == 'c':
            self.clear_lines(event)

        if event.key.lower() == 'd':
            self.delete_line(event)

        if event.key.lower() in [' ', 'enter']:
            self.accept_lines(event)

        if event.key.lower() == 'escape':
            self.refuse_spectrum(event)

    def onpick(self, event):
        """Define what happens when one clicks near the data"""
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        self.log.debug(u'Object picked at x = {}, y = {}'.format(x, y))

        i = np.argmin(np.abs(x - self.x))
        self.log.debug(u'Nearest index: {}'.format(i))

        i0 = max(0, i - 2)
        iz = min(i + 3, len(self.x))

        slice_index = np.arange(0, 5) - 2
        index = -2
        while index != 0:
            slice_index = slice_index[
                i + slice_index >= 0]  # Protect right border
            slice_index = slice_index[
                i + slice_index <= len(self.x)]  # Protect left border

            yy = self.y[i + slice_index]
            index = np.argmax(yy) - 2
            self.log.debug(u'Finding max index: {}, {}'.format(i, index))

            i = i + 1 if index > 0 else i - 1

        L = self.ax.axvline(x=self.x[i], c='red')

        self.fig.canvas.draw()
        self.vlines.append(L)
        return

    def refuse_spectrum(self, event):
        """Ignore the current spectrum and its lines"""
        plt.close('all')
        return


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
    info_fmt = u"    %(msg)s"
    wrn_fmt = TColors.WARNING + "[W]" + TColors.ENDC + \
              u" %(msg)s"
    dbg_fmt = TColors.OKBLUE + "[D]" + TColors.ENDC + \
              u" %(asctime)s %(filename)s %(funcName)s " + \
              u" %(lineno)d: %(msg)s"
    err_fmt = TColors.FAIL + "[E]" + TColors.ENDC + \
              u" %(asctime)s %(filename)s %(funcName)s " + \
              u" %(lineno)d: %(msg)s"

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
    main()