#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    Fits Combine
    
    Combine a set of FITS files using different algorithms.
    
    v0.1 - File created. No rejection here.
    
    by Quint, B. (bquint@astro.iag.usp.br)
    Dec 2013
"""

from __future__ import division, print_function

def main():
    """
    This is actually what happens when the program run.
    """
    import argparse
    import btfi
    
    # Parse arguments ---------------------------------------------------------
    parser = argparse.ArgumentParser(
             description="Combine a set of images using different algorithms.")
    parser.add_argument('-a', '--algorithm', type=str, default='average', 
                        help="Algorithm used when combining images per frame" +
                             "([average] | median | mode | sum)")
    
    parser.add_argument('-o', '--output', type=str, default="qcombined.fits", 
                        help="Name of the output cube.")
    
    parser.add_argument('-q', '--quiet', action='store_true', 
                        help="Run quietly.")
    
    parser.add_argument('-r', '--rejection', type=str, default='none', 
                        help="Type of rejection to be performed." +
                             "([none] | minmax | sigclip)")
    
    parser.add_argument('files', type=str, nargs='+', 
                        help="input filenames.")
    
    args = parser.parse_args()

    v = not args.quiet
    if v:
        print("\n Fits Combine")
        print(" by Bruno Quint (bquint@astro.iag.usp.br)")
        print(" Dec 2013 - Version 0.0.1")
        print("\n Starting program.")

    config = {'algorithm': args.algorithm,
              'lsigma': 3.0,
              'hsigma': 3.0,
              'output': args.output,
              'rejection': args.rejection,
              'verbose': not args.quiet}
    
    btfi.fits_combine(args.files, **config)

    if v:
        print("\n All done.\n")

if __name__ == '__main__':
    main()
