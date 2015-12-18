#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import division, print_function

import astropy.io.fits as pyfits
import argparse
import logging
import numpy
import os
import sys

from astropy import coordinates
from astropy import units

def main():

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(
        description="Join extensions existent in a single FITS file.")

    parser.add_argument('-a', '--acal', action='store_true',
                        help='Add astrometric calibration using the first' +
                             ' extension header as reference.')

    parser.add_argument('-b', '--bias', type=str, default=None,
                        help="Consider BIAS file for subtraction.")

    parser.add_argument('-d', '--dark', type=str, default=None,
                        help="Consider DARK file for subtraction.")

    parser.add_argument('-f', '--flat', type=str, default=None,
                        help="Consider FLAT file for division.")

    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run quietly.")

    parser.add_argument('--DEBUG', action='store_true',
                        help="Run in DEBUG mode")

    parser.add_argument('files', metavar='files', type=str, nargs='+',
                        help="input filenames.")

    args = parser.parse_args()
    v = not args.quiet

    # Create and setup log ---
    fmt = MyFormatter()
    hdlr = logging.StreamHandler(sys.stdout)
    hdlr.setFormatter(fmt)
    logging.root.addHandler(hdlr)
    logging.root.setLevel(logging.CRITICAL)
    log = logging.getLogger(__name__)

    if v:
        logging.root.setLevel(logging.INFO)

    if args.DEBUG:
        logging.root.setLevel(logging.DEBUG)

    log.info("\n SAMI - Join Extensions")
    log.info(" by Bruno Quint (bquint@astro.iag.usp.br)")
    log.info(" Mar 2015 - Version 0.4")
    log.info("\n Starting program.")

    list_of_files = args.files
    list_of_files = sorted(list_of_files)
    # number_of_files = len(list_of_files)

    log.info(" Processing data:")
    for filename in list_of_files:

        # Opening file
        prefix = "xj"
        fits_file = pyfits.open(filename)
        log.info("\n > %s " % filename)

        if args.DEBUG:
            log.debug(fits_file[0].header.__str__)
            _ = raw_input()
            log.debug(fits_file[1].header.__str__)
            _ = raw_input()
            log.debug(fits_file[2].header.__str__)
            _ = raw_input()
            log.debug(fits_file[3].header.__str__)
            _ = raw_input()
            log.debug(fits_file[4].header.__str__)
            _ = raw_input()

        # Save detector dimmensions
        w, h = str2pixels(fits_file[1].header['DETSIZE'])

        # Correct for binning
        bin_size = numpy.array(fits_file[1].header['CCDSUM'].split(' '),
                               dtype=int)
        bw, bh = w[1] // bin_size[0], h[1] // bin_size[1]

        # Create empty full frame
        new_data = numpy.empty((bh, bw), dtype=float)

        # Process each extension
        for i in range(1, 5):

            if v:
                sys.stdout.write('.')
                sys.stdout.flush()

            tx, ty = str2pixels(fits_file[i].header['TRIMSEC'])
            bx, by = str2pixels(fits_file[i].header['BIASSEC'])

            data = fits_file[i].data
            trim = data[ty[0]-1:ty[1], tx[0]-1:tx[1]]
            bias = data[by[0]-1:by[1], bx[0]-1:bx[1]]

            # Collapse the bias collumns to a single collumn.
            bias = numpy.median(bias, axis=1)

            # Fit and remove OVERSCAN
            x = numpy.arange(bias.size) + 1
            bias_fit_pars = numpy.polyfit(x, bias, 4)  # Last par = inf
            bias_fit = numpy.polyval(bias_fit_pars, x)
            bias_fit = bias_fit.reshape((bias_fit.size, 1))
            bias_fit = numpy.repeat(bias_fit, trim.shape[1], axis=1)

            trim = trim - bias_fit
            dx, dy = str2pixels(fits_file[i].header['DETSEC'])
            dx, dy = dx // bin_size[0], dy // bin_size[1]
            new_data[dy[0]:dy[1], dx[0]:dx[1]] = trim

        # Getting the main header of the FITS file instead of the header
        # an extension.
        header = fits_file[0].header

        # BIAS subtraction
        if args.bias is not None:
            bias_file = pyfits.getdata(args.bias)
            new_data = new_data - bias_file
            header['BIASFILE'] = args.bias
            prefix = 'b' + prefix

        # DARK subtraction
        if args.dark is not None:
            dark_file = pyfits.getdata(args.dark)
            new_data = new_data - dark_file
            header['DARKFILE'] = args.dark
            prefix = 'd' + prefix

        # FLAT dirvision
        if args.flat is not None:
            flat_file = pyfits.getdata(args.flat)
            new_data = new_data / flat_file
            header['FLATFILE'] = args.flat
            prefix = 'f' + prefix

        # Get astrometric calibration from the first FITS file extension.
        if args.acal is not False:

            coords = coordinates.SkyCoord(ra=fits_file[1].header['TELRA'],
                                     dec=fits_file[1].header['TELDEC'],
                                     unit=(units.hourangle, units.degree))

            header['EPOCH'] = 2000
            header['CTYPE1'] = fits_file[1].header['CTYPE1']
            header['CTYPE2'] = fits_file[1].header['CTYPE2']

            header['CRVAL1'] = fits_file[1].header['CRVAL1']
            header['CRVAL2'] = fits_file[1].header['CRVAL2']
            # header['CRVAL1'] = coords.ra.deg
            # header['CRVAL2'] = coords.dec.deg

            # header['CRPIX1'] = fits_file[1].header['CRPIX1']
            # header['CRPIX2'] = fits_file[1].header['CRPIX2']
            # header['CRPIX1'] = data.shape[0]
            # header['CRPIX2'] = data.shape[1]
            header['CRPIX1'] = 0
            header['CRPIX2'] = 0

            header['CD1_1'] = fits_file[1].header['CD1_1'] # * bin_size[0]
            header['CD2_1'] = fits_file[1].header['CD2_1'] # * bin_size[0]
            header['CD1_2'] = fits_file[1].header['CD1_2'] # * bin_size[0]
            header['CD2_2'] = fits_file[1].header['CD2_2'] # * bin_size[0]
            header['CDELT1'] = fits_file[1].header['CDELT1'] # * bin_size[0]
            header['CDELT2'] = fits_file[1].header['CDELT2'] # * bin_size[0]

            prefix = 'a' + prefix

        # Removing bad column and line
        n_rows, n_columns = new_data.shape
        temp_column = new_data[:, n_columns//2-1:n_columns//2+1]
        new_data[:, n_columns//2-1:-2] = new_data[:, n_columns//2+1:]
        new_data[:, -2:] = temp_column

        # Writing file
        path, filename = os.path.split(filename)
        pyfits.writeto(os.path.join(path, prefix + filename), new_data, header,
                       clobber=True)

    log.info("\n All done!")


def str2pixels(my_string):

    my_string = my_string.replace('[', '')
    my_string = my_string.replace(']', '')
    x, y = my_string.split(',')

    x = x.split(':')
    y = y.split(':')

    # "-1" fix from IDL to Python
    x = numpy.array(x, dtype=int)
    y = numpy.array(y, dtype=int)

    return x, y


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
