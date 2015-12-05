#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    Roll data-cube in Z.

    This script was created to easily roll a data-cube in Z keeping the
    wavelength calibration.
"""
from __future__ import print_function, division
from astropy.io import fits as pyfits
from numpy import roll
from os.path import join, split, splitext

import argparse


__author__ = 'b1quint'

def main():

    ## Parsing arguments ---
    args = parse_arguments()
    v = not args.quiet

    ## Print the script's header?
    print_header(verbose=v)

    ## Load data ---
    data, header = load_data(args.filename, verbose=v)

    ## Roll data --
    data = roll(data, args.shift_size, axis=0)
    if args.keep_calibration:
        try:
            header['CRPIX3'] = header['CRPIX3'] + args.shift_size
        except KeyError:
            header['CRPIX3'] = 1 + args.shift_size

    ## Prepare output filename --
    if args.output is None:
        path, output = split(args.filename)
        output, extension = splitext(output)
        output = join(path, output + '.roll' + extension)
    else:
        output = args.output

    output = safe_save(output, verbose=v, overwrite=True)

    ## Write output data --
    pyfits.writeto(output, data, header, clobber=True)

    return


def load_data(filename, verbose=False):
    """
    Load the data and the header of the input file.

    :param filename: a string containing the path for the input data-cube.
    :arg verbose: verbose mode (true/False)?
    :return data: a numpy.ndarray containing the 3D data-cube.
    :return header: a pyfits.Header instance containing the header of the file.
    """
    v = verbose
    if v: print(" Loading data from file: %s" % filename)
    data = pyfits.getdata(filename)
    header = pyfits.getheader(filename)
    if v: print(" Done!")

    return data, header


def parse_arguments():
    """
    Parse the arguments given by the user.
    :return:
    """
    parser = argparse.ArgumentParser(description="Fits an existing phase-map.")

    parser.add_argument('filename',
                        type=str,
                        help="Input data-cube name.")

    parser.add_argument('shift_size',
                        type=int,
                        help="Shift size in number of channels.")

    parser.add_argument('-k', '--keep_calibration', action='store_true',
                        help="Keep the calibration while rolling?")

    parser.add_argument('-o', '--output', type=str, default=None,
                        help="Name of the output data-cube.")

    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run program quietly.")

    return parser.parse_args()


def print_header(verbose=False):
    """
    Simply prints a header if verbose mode is true.
    :param verbose:
    :return: None.
    """
    if verbose:
        print("")
        print(" Roll Data-Cube")
        print(" by Bruno C. Quint - Apr 2015")
        print(" Version 0.0a")
        print("")

    return

def safe_save(name, overwrite=False, verbose=False):
    """
    This is a generic method used to check if a file called 'name' already
    exists. If so, it starts some interaction with the user.

    @param name: the name of the file that will be written in the future.

    @keyword overwrite: if False, this method will interact with the user to
    ask if 'name' file shall be overwritten or if a new name will be given. If
    True, 'name' file is automatically overwritten.

    @keyword verbose: force verbose mode on even when overwrite is automatic.

    v1.0.1 - added 'overwrite' keyword.
           - added 'verbose' keyword.
    """
    import os
    import sys

    v = False if (overwrite is True) else True
    if v:
        print("\n Writing to output file %s" % name)

    while os.path.exists(name):

        if overwrite in ['y', 'Y', True]:
            if v or verbose:
                print(" Overwriting %s file.\n" % name)
            os.remove(name)

        elif overwrite in ['', 'n', 'N', False]:
            name = raw_input("   Please, enter a new filename:\n   > ")

        elif overwrite in ['q']:
            if v:
                print(" Exiting program.\n")
            sys.exit()

        else:
            overwrite = raw_input("   '%s' file exist. Overwrite? (y/[n])"% name)
            if v:
                print(" Writing data-cube to %s\n" % name)

    return name

if __name__ == '__main__':
    main()