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
import os
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
        description=u'Calculate the wavelength calibration for FP observation.'+
                    u' The parameters not given in the command line will be '+
                    u'asked during runtime. Wavelength 1 has to have two ' +
                    u'identified peaks and Wavelength 2 needs one identified '+
                    u'peak.')

    parser.add_argument('-v', '--verbose_level', default=1, type=int,
                        help=u'Set verbose level: 0 - Quiet, 1 - Info,' +
                             u"2 - Warnings, 3 - Debug. (Default: 1)")
    parser.add_argument('--w1', default=None, type=float,
                        help=u"First wavelength in angstroms.")
    parser.add_argument('--w2', default=None, type=float,
                        help=u"Second wavelength in angstroms.")
    parser.add_argument("--file1", default=None, type=str,
                        help=u"Input FITS file for calibration containing w1.")
    parser.add_argument("--file2", default=None, type=str,
                        help=u"Input FITS file for calibration containing w2.")

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

    # Gather missing information ---
    log.info(u"\n    Fabry-Perot Wavelength Calibration\n" + \
             u"    by Bruno C. Quint - bquint at ctio dot noao dot edu\n" + \
             u"    Version 0.0a\n")

    log.info(u"For the wavelength calibration you will need two wavelengths.\n"+
             u"    One wavelength has to have two identified peaks.\n" +
             u"    The other will have to have one identified peak.\n" +
             u"    Press <ENTER> when you are ready.\n")
    dummy = None
    while dummy is None:
        dummy = raw_input()

    if args.w1 == None:
        log.info(u"Now, please enter the first known wavelength in [A]:")
        wavelength_1 = raw_input("  > ")
    else:
        wavelength_1 = args.w1

    if args.w2 == None:
        log.info(u"Now, enter the second known wavelength in [A]:")
        wavelength_2 = raw_input("  > ")
    else:
        wavelength_2 = args.w2

    if args.file1 == None:
        log.info(u"Enter the path to the FITS file that contains a spectrum" +
                 u" with the first wavelength.")
        filename_1 = raw_input(" > ")
    else:
        filename_1 = args.file1

    if args.file2 == None:
        log.info(u"Enter the path to the FITS file that contains a spectrum" +
                 u" with the second wavelength.")
        filename_2 = raw_input(" > ")
    else:
        filename_2 = args.file2

    # Make sure that the wavelengths are float
    wavelength_1 = float(wavelength_1)
    wavelength_2 = float(wavelength_2)

    # Process first file ---
    f1 = None
    while f1 == None:
        log.info(u"Reading file {}".format(filename_1))
        try:
            f1 = pyfits.open(filename_1)[0]
        except IOError:
            log.error(u"{} not found.".format(filename_1))
            log.error(u"Enter a new filename: ")
            filename_1 = raw_input(" > ")
            continue

    try:
        z1 = (np.arange(f1.header['NAXIS1']) - f1.header['CRPIX1'] + 1) * \
            f1.header['CDELT1'] + f1.header['CRVAL1']
    except KeyError:
        log.warning(u"Spectral calibration not found in " +
                    u"{}'s header.".format(filename_1))
        z1 = np.arange(f1.header['NAXIS1'])
        f1.header['CRPIX1'] = 1
        f1.header['CRVAL1'] = 1
        f1.header['CDELT1'] = 1

    s1 = f1.data
    s1 = np.convolve(s1, [0.1, 0.15, 0.5, 0.15, 0.1], mode='same')

    tck1 = splrep(z1, s1, k=5)
    znew1 = np.linspace(z1[0], z1[-1], num=1000, endpoint=True)
    snew1 = splev(znew1, tck1, der=0)

    # Find first peak for the first wavelength ---
    log.info(u"\n    A new window has opened. Please, click near the peak" + \
             u" that you want to identify (first peak for first wavelength).")
    log.info(u"Click near a peak to add it.")
    log.info(u"Press:")
    log.info(u"<space_bar> or <ENTER> to accept the line selected.")
    log.info(u'<Q> to leave the program.')

    p = SingleLinePicker(filename_1, z1, s1, znew1, snew1)
    z11 = p.xs
    log.info(u"First peak found at {0:0.2f}".format(z11) +
             u" for wavelength_1 = {0:0.2f} A".format(wavelength_1))

    # Find second peak for the first wavelength ---
    log.info(u"\n    A new window has opened. Please, click near the peak" + \
             u" that you want to identify (more to the right).")
    log.info(u"Click near a peak to add it.")
    log.info(u"Press:")
    log.info(u"<space_bar> or <ENTER> to accept the line selected.")
    log.info(u'<Q> to leave the program.')

    p = SingleLinePicker(filename_1, z1, s1, znew1, snew1)
    z12 = p.xs
    log.info(u"Second peak found at {0:0.2f}".format(z12) +
             u" for wavelength_1 = {0:0.2f} A".format(wavelength_1))

    # Process second file ---
    f2 = None
    while f2 == None:
        log.info(u"Reading file {}".format(filename_2))
        try:
            f2 = pyfits.open(filename_2)[0]
        except IOError:
            log.error(u"{} not found.".format(filename_2))
            log.error(u"Enter a new filename: ")
            filename_1 = raw_input(" > ")
            continue

    try:
        z2 = (np.arange(f2.header['NAXIS1']) - f2.header['CRPIX1'] + 1) * \
            f2.header['CDELT1'] + f2.header['CRVAL1']
    except KeyError:
        log.warning(u"Spectral calibration not found in " +
                    u"{}'s header.".format(filename_1))
        z2 = np.arange(f1.header['NAXIS1'])
        f2.header['CRPIX1'] = 1
        f2.header['CRVAL1'] = 1
        f2.header['CDELT1'] = 1

    s2 = f2.data
    s2 = np.convolve(s2, [0.1, 0.15, 0.5, 0.15, 0.1], mode='same')

    tck2 = splrep(z2, s2, k=5)
    znew2 = np.linspace(z2[0], z2[-1], num=1000, endpoint=True)
    snew2 = splev(znew2, tck2, der=0)

    log.info(u"\n    A new window has opened. Please, click near the peak" + \
             u" that you want to identify (more to the left).")
    log.info(u"Click near a peak to add it.")
    log.info(u"Press:")
    log.info(u"<space_bar> or <ENTER> to accept the line selected.")
    log.info(u'<Q> to leave the program.')

    p = SingleLinePicker(filename_2, z2, s2, znew2, snew2)
    z21 = p.xs
    log.info(u"Peak found at {0:0.2f}".format(z21) +
             u" for wavelength_2 = {0:0.2f} A".format(wavelength_2))

    # Calculating the wavelength calibration ---
    curvature = 0
    while curvature not in [1, 2]:
        curvature = raw_input(u"\n    Do the rings increase or decrease with Z?" +
                              u"\n    0 - Increase\n    1 - Decrease\n > ")
        try:
            curvature = float(curvature)
        except ValueError:
            pass

    dz1 = np.abs(z11 - z12)
    B = ((-1) ** (curvature) * wavelength_1 / dz1)
    log.info("B = {0:+0.5f} A / bnv".format(B))

    d = 200. # [um]
    d *= 1e4  # [um] to [A]
    m1 = 2 * d / wavelength_1
    m1 = int(m1) + np.arange(50) - 25

    c = curvature
    m2 = (wavelength_1 / wavelength_2) * (m1 - (-1) ** c * (z11 - z21) / dz1)
    distance_from_integer = np.abs(m2 - np.round(m2, 0))

    i = np.argmin(distance_from_integer)
    m1 = m1[i]
    m2 = m2[i]

    A = m1 * wavelength_1 - B * z11
    log.info("A = {0:+0.5f} A".format(A))

    CRPIX1 = f1.header['CRPIX1']
    CRVAL1 = (A + B * f1.header['CRVAL1']) / m1
    CDELT1 = (B * f1.header['CDELT1']) / m1

    log.info(u"\n    Calibration parameters found for" +
             u" {0:0.1f}A".format(wavelength_1) +
             u" at order {0:0d}".format(m1))
    log.info(u"CRPIX1 = {}".format(CRPIX1))
    log.info(u"CRVAL1 = {}".format(CRVAL1))
    log.info(u"CDELT1 = {}".format(CDELT1))


    w = (np.arange(z1.size) - CRPIX1 + 1) * \
            CDELT1 + CRVAL1
    log.info(u"For this configuration the wavelength range is:")
    log.info(u"{0:0.2f} A - {1:0.2f} A".format(w[0], w[-1]))
    log.info(u"Type one of the following options and press enter for " +
             u"different actions:" +
             u"\n    '+' If you want to increase an order press" +
             u"\n    '-' If you want to decrease an order press" +
             # u"\n    'a' If you want to apply the current calibration to " +
             # u"a FITS file." +
             u"\n    'q' to leave" +
             u"\n    A channel number, if you want to evaluate the " +
             u"wavelength at that channel."
             '')

    a = None
    while a not in ['q']:
        a = raw_input(' > ')

        if a == '+':
            m1 += 1

        elif a == '-':
            m1 -= 1

        elif a.isdigit():

            w = apply_calibration_to_channel(float(a),
                {'CRPIX1': f1.header['CRPIX1'], 'CRVAL1': f1.header['CRVAL1'],
                 'CDELT1': f1.header['CDELT1']},
                wcal)

            log.info(u"Wavelength at {} = {:0.2f} A\n".format(a, w))
            _ = raw_input(u"    Press <ENTER> to continue.")

        CRPIX1 = f1.header['CRPIX1']
        CRVAL1 = (A + B * f1.header['CRVAL1']) / m1
        CDELT1 = (B * f1.header['CDELT1']) / m1
        wcal = {'CRPIX1': CRPIX1, 'CRVAL1': CRVAL1,'CDELT1': CDELT1}

        log.info(u"\n    Calibration parameters found for" +
                 u" {0:0.1f}A".format(wavelength_1) +
                 u" at order {0:0d}".format(m1))
        log.info(u"CRPIX1 = {}".format(CRPIX1))
        log.info(u"CRVAL1 = {}".format(CRVAL1))
        log.info(u"CDELT1 = {}".format(CDELT1))


        w = (np.arange(z1.size) - CRPIX1 + 1) * \
                CDELT1 + CRVAL1
        log.info(u"For this configuration the wavelength range is:")
        log.info(u"{0:0.2f} A - {1:0.2f} A".format(w[0], w[-1]))
        log.info(u"Reference wavelength is: {0:.2f}".format(wavelength_1))
        log.info(u"Reference order is: {0:d}".format(m1))
        log.info(u"Type one of the following options and press enter for " +
                 u"different actions:" +
                 u"\n    '+' If you want to increase an order press" +
                 u"\n    '-' If you want to decrease an order press" +
                 u"\n    'q' to leave" +
                 u"\n    A channel number, if you want to evaluate the " +
                 u"wavelength at that channel."
                 '')

    print('')
    return 0


def apply_calibration_to_channel(z, old_cal, new_cal):

    channel = (z - old_cal['CRVAL1']) / old_cal['CDELT1'] + old_cal['CRPIX1'] -1

    wavelength = (channel - new_cal['CRPIX1'] + 1) * new_cal['CDELT1'] + \
        new_cal['CRVAL1']

    return wavelength


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

    def accept_lines(self):
        """Accept the current lines drawn"""
        for line in self.ax.lines[2:]:
            x, _ = line.get_xdata()
            self.xs.append(x)
        plt.close('all')
        return

    def clear_lines(self):
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
            self.clear_lines()

        if event.key.lower() == 'd':
            self.delete_line(event)

        if event.key.lower() in [' ', 'enter']:
            self.accept_lines()

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


class SingleLinePicker(object):

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

    def accept_lines(self):
        """Accept the current lines drawn"""

        self.xs = np.average(self.ax.lines[2].get_xdata())
        plt.close('all')

        return

    def delete_line(self, event):
        """Delete a single line near the cursor."""
        for line in self.ax.lines[2:]:
            xdata, _ = line.get_xdata()
            self.log.debug(u"xdata = {0:0.2f}".format(xdata) +
                           u", event.xdata" +
                           " = {0:0.2f}".format(event.mouseevent.xdata))
            xdata_index = np.argmin(np.abs(self.x - xdata))
            exdata_index = np.argmin(np.abs(self.x - event.mouseevent.xdata))
            self.log.debug(u"xdata_i = {:d}".format(xdata_index) +
                           u", event.xdata_i" +
                           u" = {:d}".format(exdata_index))
            try:
                line.remove()
            except ValueError:
                pass
            self.fig.canvas.draw()
        return

    def onpick(self, event):
        """Define what happens when one clicks near the data"""

        self.delete_line(event)

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

    def onpress(self, event):
        """define some key press events"""
        self.log.debug(u"Key pressed: " + event.key.lower())

        if event.key == 'Q':
            self.log.info(u'You pressed "Q". Leaving the program now.')
            self.log.info(u'Bye!')
            sys.exit()

        if event.key.lower() in [' ', 'enter']:
            self.accept_lines()

        return


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


# Custom formatter
class MyFormatter(logging.Formatter):
    info_fmt = u"    %(msg)s"
    wrn_fmt = TColors.WARNING + "[W]" + TColors.ENDC + \
        u" %(msg)s"
    dbg_fmt = u"{0}[D]{1} %(asctime)s %(filename)s %(funcName)s  %(lineno)d: %(msg)s".format(
        TColors.OKBLUE, TColors.ENDC)
    err_fmt = u"{0}[E]{1} %(asctime)s %(filename)s %(funcName)s  %(lineno)d: %(msg)s".format(
        TColors.FAIL, TColors.ENDC)

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
