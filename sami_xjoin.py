#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import division, print_function

import astropy.io.fits as pyfits
import argparse
import numpy as np
import os
import sys

bad_columns = [
    [167, 0, 512],
    [476, 0, 513],
    [602, 0, 513],
    [671, 0, 513],
    [810, 0, 513],
    [918, 0, 513],
    [917, 0, 513],
    [213, 513, 1024]
    ]

bad_lines = [
    [214, 239, 688],
    [477, 516, 490],
    [387, 429, 455],
    [574, 603, 494],
    [574, 603, 493],
    [640, 672, 388],
    [604, 671, 388]
    ]


def main():

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(
            description="Join extensions existent in a single FITS file."
    )

    parser.add_argument('-b','--bias', type=str, default=None,
                        help="Consider BIAS file for subtraction.")
    parser.add_argument('-c','--clean', action='store_true',
                        help="Clean known bad columns and lines by taking the "
                             "median value of their neighbours.")
    parser.add_argument('-d','--dark', type=str, default=None,
                        help="Consider DARK file for subtraction.")
    parser.add_argument('-f','--flat', type=str, default=None,
                        help="Consider FLAT file for division.")
    parser.add_argument('-g','--glow', type=str, default=None,
                        help="Consider DARK file to correct lateral glows.")
    parser.add_argument('-t','--exptime', action='store_true',
                        help="Divide by exposure time.")
    parser.add_argument('-q','--quiet', action='store_true',
                        help="Run quietly.")
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                        help="input filenames.")

    args = parser.parse_args()
    v = not args.quiet
    if v:
        print("\n SAMI - Join Extensions")
        print(" by Bruno Quint (bquint@astro.iag.usp.br)")
        print(" Mar 2015 - Version 0.4")
        print("\n Starting program.")

    list_of_files = args.files
    list_of_files = sorted(list_of_files)
    number_of_files = len(list_of_files)

    if v: print(" Processing data:")
    for filename in list_of_files:

        if v:
            sys.stdout.write("\n > %s " % filename)
            sys.stdout.flush()

        prefix = "xj"
        fits_file = pyfits.open(filename)
        w, h = str2pixels(fits_file[1].header['DETSIZE'])

        # Correct for binning
        bin_size = np.array(fits_file[1].header['CCDSUM'].split(' '), dtype=int)
        bw, bh = w[1] // bin_size[0], h[1] // bin_size[1]

        # Create empty full frame
        new_data = np.empty((bh,bw), dtype=float)

        # Process each extension
        for i in range(1, 5):

            if v:
                sys.stdout.write('.')
                sys.stdout.flush()

            tx, ty = str2pixels(fits_file[i].header['TRIMSEC'])
            bx, by = str2pixels(fits_file[i].header['BIASSEC'])

            data = fits_file[i].data
            trim = data[ty[0]-1:ty[1],tx[0]-1:tx[1]]
            bias = data[by[0]-1:by[1],bx[0]-1:bx[1]]

            # Collapse the bias columns to a single column.
            bias = np.median(bias, axis=1)

            # Fit and remove OVERSCAN
            x = np.arange(bias.size) + 1
            bias_fit_pars = np.polyfit(x, bias, 4) # Last par = inf
            bias_fit = np.polyval(bias_fit_pars, x)
            bias_fit = bias_fit.reshape((bias_fit.size, 1))
            bias_fit = np.repeat(bias_fit, trim.shape[1], axis=1)

            trim = trim - bias_fit
            dx, dy = str2pixels(fits_file[i].header['DETSEC'])
            dx, dy = dx // bin_size[0], dy // bin_size[1]
            new_data[dy[0]:dy[1],dx[0]:dx[1]] = trim

        # Getting the main header of the FITS file instead of the header
        # an extension.
        header = fits_file[0].header
        header.append('UNITS')
        header.set('UNITS', value='COUNTS', comment='Pixel intensity units.')

        # BIAS subtraction
        if args.bias is not None:
            bias_file = pyfits.getdata(args.bias)
            new_data = new_data - bias_file
            header['BIASFILE'] = args.bias
            header.add_history('Bias subtracted')
            prefix = 'b' + prefix

        # DARK subtraction
        if args.dark is not None:
            dark_file = pyfits.getdata(args.dark)
            new_data = new_data - dark_file
            header['DARKFILE'] = args.dark
            prefix = 'd' + prefix
            header.add_history('Dark subtracted')

        # FLAT division
        if args.flat is not None:
            flat_file = pyfits.getdata(args.flat)
            new_data = new_data / flat_file
            header['FLATFILE'] = args.flat
            header.add_history('Flat normalized') 
            prefix = 'f' + prefix

        # Normalize by the EXPOSURE TIME
        if args.exptime is True:
            try:
                exptime = float(header['EXPTIME'])
                new_data /= exptime
                header['UNITS'] = 'COUNTS/s'
                header.add_history('Divided by exposure time.')
                prefix = 't' + prefix
            except KeyError:
                pass

        # Removing bad column and line
        n_rows, n_columns = new_data.shape
        temp_column = new_data[:,n_columns//2-1:n_columns//2+1]
        new_data[:,n_columns//2-1:-2] = new_data[:,n_columns//2+1:]
        new_data[:,-2:] = temp_column

        # Clean known bad columns and lines
        if args.clean is True:
            new_data = clean_columns(new_data)
            # new_data = clean_lines(new_data)
            header.add_history('Cleaned bad columns and lines.')
            prefix = 'c' + prefix

        # Remove lateral glows
        if args.glow is not None:
            dark = pyfits.getdata(args.glow)
            new_data = remove_glows(new_data, dark)
            header.add_history('Lateral glow removed using %s file' % args.dark)
            prefix = 'g' + prefix

        # Writing file
        header.add_history('Extensions joined using "sami_xjoin"')
        path, filename = os.path.split(filename)
        pyfits.writeto(os.path.join(path, prefix + filename), new_data, header, clobber=True)

    print("\n All done!")


def clean_column(_data, x0, y0, yf, n=5):
    t1 = _data[y0:yf, x0 - n:x0]
    t2 = _data[y0:yf, x0 + 1:x0 + n]
    t = np.hstack((t1, t2))
    _data[y0:yf, x0] = np.median(t)
    return _data


def clean_columns(_data):
    global bad_columns
    for column in bad_columns:
        x0 = column[0]
        y0 = column[1]
        yf = column[2]
        _data = clean_column(_data, x0, y0, yf)
    return _data


def clean_line(_data, x0, xf, y, n=5):
    t1 = _data[y - n:y, x0:xf]
    t2 = _data[y + 1:y + 1, x0:xf]
    t = np.vstack((t1, t2))
    _data[y, x0:xf] = np.median(t)
    return _data


def clean_lines(_data):
    global bad_lines
    for line in bad_lines:
        x0 = line[0]
        xf = line[1]
        y = line[2]
        _data = clean_line(_data, x0, x0, y)
    return _data


def remove_glows(_data, _dark):

    _dark = clean_columns(_dark)
    _dark = clean_lines(_dark)

    dark_midpt1 = np.median(_dark[539:589, 999:1009])
    dark_midpt2 = np.median(_dark[449:506, 975:1019])
    dark_diff = dark_midpt2 - dark_midpt1
    _dark -= dark_midpt1

    midpt1 = np.median(_data[539:589, 999:1009])
    midpt2 = np.median(_data[449:506, 975:1019])
    diff = midpt2 - midpt1

    k = diff / dark_diff
    temp_dark = _dark * k
    _data -= midpt1
    _data -= temp_dark

    return _data


def str2pixels(my_string):

    my_string = my_string.replace('[','')
    my_string = my_string.replace(']','')
    x, y = my_string.split(',')

    x = x.split(':')
    y = y.split(':')

    # "-1" fix from IDL to Python
    x = np.array(x, dtype=int)
    y = np.array(y, dtype=int)

    return x, y

if __name__ == '__main__':
    main()
