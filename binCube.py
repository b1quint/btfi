#!/usr/bin/python
#!-*- coding: utf8 -*-

from __future__ import division, print_function;

def binCube(input_cubename, bin_size, output_cubename=None, verbose=False, keep_header=False):
    """
    Apply spatial binning on a data-cube.
    
    @param input_cubename:    path containing the name of the input FITS data-cube.
    @param bin_size:          bin size in pixel units.
    @keyword keep_header:     Keep header untouched?
    @keyword output_cubename: path containing the name of the output binned FITS data-cube.
    @keyword verbose:         turn on verbose mode.
    """
    import numpy;
    import os;
    import pyfits;
    
    v = verbose;
    if v: print(" Loading the following file:\n  %s" % input_cubename);
    cube = pyfits.getdata(input_cubename);
    
    if v: print(" Allocating memory.");
    temp = numpy.zeros((cube.shape[0], cube.shape[1] // bin_size, cube.shape[2] // bin_size));
    
    if v: print(" Binning cube.");
    for i in range(bin_size):
        temp += cube[:,i::bin_size,i::bin_size];
    temp = temp / bin_size;
    
    if not keep_header:
        if v: print(" Fixing header.");
        header = pyfits.getheader(input_cubename);
        try: 
            header['CDELT1'] = header['CDELT1'] * bin_size;
            header['CDELT2'] = header['CDELT2'] * bin_size;
        except KeyError:
            header['CDELT1'] = bin_size;
            header['CDELT2'] = bin_size;

    if v: print(" Writing output file");        
    if output_cubename is None:
        output_cubename = os.path.splitext(input_cubename)[0] + "_binned.fits";
        output_cubename = safesave(output_cubename);
    pyfits.writeto(output_cubename, temp, header);
         
def safesave(name, overwrite=False, verbose=False):
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
    import os;
    import sys;

    v = False if (overwrite is True) else True;
    if v: print("  Writing to output file %s" % name);

    while os.path.exists(name):

        if overwrite in ['y', 'Y', True]:
            if v or verbose: print("   Overwriting %s file." % name);
            os.remove(name);

        elif overwrite in ['', 'n', 'N', False]:
            name = raw_input("   Please, enter a new filename:\n   > ");

        elif overwrite in ['q']:
            if v: print("   Exiting program.");
            sys.exit();

        else:
            overwrite = raw_input("   '%s' file exist. Overwrite? (y/[n])"%name);
            if v: print("   Writing data-cube to %s" %name);

    return name        
        

if __name__ == '__main__':
    
    import argparse;
    import time;
    
    parser = argparse.ArgumentParser(description="Apply spatial binning on a data-cube using the average of the pixels inside the bin window.");
    parser.add_argument('-o','--output', type=str, default=None, help="Name of the output binned data-cube.");
    parser.add_argument('-k','--keep', action='store_true', help="Keep header untouched?");
    parser.add_argument('-q','--quiet', action='store_true', help="Run it quietly.");
    parser.add_argument('input_cube', metavar='input_cube', type=str, help="Input data-cube to be binned.");
    parser.add_argument('bin_size', metavar='bin_size', type=int, help="Bin size.");
    args = parser.parse_args();
    v = not args.quiet;
    
    if v: start = time.time();
    if v: print("\n bin_datacube.py - by Bruno C. Quint");
    if v: print(" September 2013 - v1.0");
    
    binCube(args.input_cube, args.bin_size, output_cubename=args.output, verbose=v, keep_header=args.keep);
    
    if v: end = time.time() - start;
    if v: print("\n Total time ellapsed: %02d:%02d:%02d" % (end // 3600, end % 3600 // 60, end % 60));     
    if v: print(" All done!\n");
