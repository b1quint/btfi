#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    Data-Cube - Change Z Calibration
"""

from __future__ import division, print_function

import astropy.io.fits as pyfits
import argparse

__author__ = 'Bruno Quint'
__date__ = '2015.12.16'
__ver__ = '0.0.0'


class Main(object):
    def __init__(self):
        return

    def run(self, options):

        # Open file
        print(options.filename)
        data = pyfits.getdata(options.filename)
        header = pyfits.getheader(options.filename)

        # Check information content
        if options.reference is None:
            options.reference = ''
            while not options.reference.isdigit():
                options.reference = raw_input(" Please, enter reference channel:\n ")
                
        if options.value is None:
            options.value = ''
            while not options.value.isdigit():
                options.value = raw_input(" Please, enter value channel:\n ")
                
        if options.delta is None:
            options.delta = ''
            while not options.delta.isdigit():
                options.delta = raw_input(" Please, enter delta channel:\n ")
                
        if options.units is None:
            options.units = '0.0'
            while not options.units.isalpha():
                options.units = raw_input(" Please, enter units channel:\n ")

        # Update header
        header['CRPIX3'] = int(options.reference)
        header['CRVAL3'] = float(options.value)
        header['CDELT3'] = float(options.delta)
        header['C3_3'] = float(options.delta)
        header['CUNIT3'] = options.units

        # Write file
        options.filename = self.safe_save(options.filename, overwrite=True)
        pyfits.writeto(options.filename, data, header)

        return

    def safe_save(self, name, overwrite=False, verbose=False):
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
            print("\n Writing to output file %s" % name)

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
                overwrite = raw_input("   '%s' file exist. Overwrite? (y/[n])" %
                                      name)
                if v:
                    print(" Writing data-cube to %s" % name)

        return name

if __name__ == '__main__':

    # Parsing Arguments
    parser = argparse.ArgumentParser(
        description="This script adds/changes the calibration of the third" +
                    " axis of a data-cube.")

    parser.add_argument('filename', type=str, help="input filename.")

    parser.add_argument('-r', '--reference', type=float, default=None,
                        help="Set reference channel.")

    parser.add_argument('-a', '--value', type=float, default=None,
                        help="Set reference value.")

    parser.add_argument('-d', '--delta', type=float, default=None,
                        help="Set increment between channels.")

    parser.add_argument('-u', '--units', type=str, default=None,
                        help="Set units")

    parser.add_argument('-v','--verbose', type=int, default=0,
                        help="Verbose level: 0 - None, 1 - Info, 2 - Debug.")

    args = parser.parse_args()

    main = Main()
    main.run(args)
