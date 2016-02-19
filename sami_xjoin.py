#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import division, print_function

import astropy.io.fits as pyfits
import argparse
import numpy
import os
import sys

def main():

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(description="Join extensions existent in a single FITS file.")
    parser.add_argument('-b','--bias', type=str, default=None, help="Consider BIAS file for subtraction.")
    parser.add_argument('-d','--dark', type=str, default=None, help="Consider DARK file for subtraction.")
    parser.add_argument('-f','--flat', type=str, default=None, help="Consider FLAT file for division.")
    parser.add_argument('-q','--quiet', action='store_true', help="Run quietly.")
    parser.add_argument('files', metavar='files', type=str, nargs='+', help="input filenames.")

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
        bin_size = numpy.array(fits_file[1].header['CCDSUM'].split(' '), dtype=int)
        bw, bh = w[1] // bin_size[0], h[1] // bin_size[1]

        # Create empty full frame
        new_data = numpy.empty((bh,bw), dtype=float)

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
            bias = numpy.median(bias, axis=1) # Collapse the bias collumns to a single collumn.

            # Fit and remove OVERSCAN
            x = numpy.arange(bias.size) + 1
            bias_fit_pars = numpy.polyfit(x, bias, 4) # Last par = inf
            bias_fit = numpy.polyval(bias_fit_pars, x)
            bias_fit = bias_fit.reshape((bias_fit.size, 1))
            bias_fit = numpy.repeat(bias_fit, trim.shape[1], axis=1)

            trim = trim - bias_fit
            dx, dy = str2pixels(fits_file[i].header['DETSEC'])
            dx, dy = dx // bin_size[0], dy // bin_size[1]
            new_data[dy[0]:dy[1],dx[0]:dx[1]] = trim

        # Getting the main header of the FITS file instead of the header
        # an extention.
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

        # Removing bad column and line
        n_rows, n_columns = new_data.shape
        temp_column = new_data[:,n_columns//2-1:n_columns//2+1]
        new_data[:,n_columns//2-1:-2] = new_data[:,n_columns//2+1:]
        new_data[:,-2:] = temp_column

        # Writing file
        path, filename = os.path.split(filename)
        pyfits.writeto(os.path.join(path, prefix + filename), new_data, header, clobber=True)

    print("\n All done!")

def str2pixels(my_string):

    my_string = my_string.replace('[','')
    my_string = my_string.replace(']','')
    x, y = my_string.split(',')

    x = x.split(':')
    y = y.split(':')

    # "-1" fix from IDL to Python
    x = numpy.array(x, dtype=int)
    y = numpy.array(y, dtype=int)

    return x, y

if __name__ == '__main__':
    main()
