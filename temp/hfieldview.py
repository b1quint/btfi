#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import astropy.io.fits as pyfits

def main():

    # Setting Options ---------------------------------------------------------
    parser = argparse.ArgumentParser(description="Read and print a fits" + \
                                     " file header (full or field)" )
                        
    parser.add_argument('field', type=str, help="Field to be printed")
    
    parser.add_argument('value', type=str, help="Value to be updated")
    
    parser.add_argument('fits_files', type=str, nargs='+', 
                        help="Input fits files.")
    
    args = parser.parse_args()
    print("\n# Header Field View - by Bruno Quint - v0.0")
    
    # Reading header ----------------------------------------------------------
    print('# %-45s \t Field="%s"' % ("Filename", args.field))
    for fits_file in args.fits_files:
        header = pyfits.getheader(fits_file)
        try:
            value = header[args.field]
        except(KeyError):
            value = '-'
        print("  %-45s \t %s" % (fits_file, value))
    print("")
        
if __name__ == '__main__':
    main()
