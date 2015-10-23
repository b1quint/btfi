#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    hcopy.py
    
    This file contains methods used to copy one FITS file's header to another(s).
    
    by Bruno C. Quint (bquint@astro.iag.usp.br)
    August 2013
"""
from __future__ import division, print_function


def hcopy(input_file, output_files, **kwargs):
    """
    Copy the header or part of it from one file to another.
    
    @param input_file: Input file containing the original header.
    @param output_files: File(s) that will receive the new header/field.
    
    @keyword field: Name of the field that will be copied (in the case of a particular field).
    @keyword quiet: Turn verbose mode off. 
    """
    from os.path import exists
    from astropy.io.fits import getheader, getdata, writeto
    
    import warnings
    warnings.filterwarnings('ignore')
    
    verbose = True if 'quiet' not in kwargs else not kwargs['quiet'] 
    
    if not exists(input_file): raise ValueError, "Input File does not exist.";
    
    if output_files.__class__ is ''.__class__: output_files = [output_files];
    elif output_files.__class__ is [''].__class__: pass;
    else: raise TypeError, "'output_file' must be a string or a list of strings.";
    
    if verbose: print("\n Reading header from %s file." % input_file);
    header = getheader(input_file);
    if verbose: print(" Done.")
    
    for output_file in output_files:
       
        if verbose: print("\n Reading %s file." % output_file);    
        odata = getdata(output_file);
        oheader = getheader(output_file);
        
        if 'field' in kwargs and kwargs['field'] is not None:
            field = kwargs['field'];
            if verbose: print(" Updating %s field." % field);
            oheader.set(field, header[field], header.comments[field])
        else:
            if verbose: print(" Updating header.");
            oheader = header;
        
        if verbose: print (" Saving file.");
        writeto(output_file, odata, oheader, clobber=True)
        
        if verbose: print(" Done.");
    
    return

if __name__ == '__main__':

    import argparse
    import time
    
    t = time.time()    
    parser = argparse.ArgumentParser(
        description="Copy the header or a header's field from one FITS file " +
                    "to other(s)")
    parser.add_argument('--field', '-f', type=str, default=None,
                        help="The particular field that will be copied.")
    parser.add_argument('--quiet', '-q', action='store_true',
                        help="Turn off verbose mode")
    parser.add_argument('input', type=str,
                        help="Input file containing the original header.")
    parser.add_argument('output', type=str, nargs='+',
                        help="Files that will receive the new header/field")
    args = parser.parse_args()
    
    if not args.quiet: 
        print("\n HCOPY - " +
              "Copy information from one FITS file's header to another.")
        print(" by Bruno C. Quint - bquint at astro.iag.usp.br")
        print(" v0.1 2013.08")
    
    hcopy(args.input, args.output, field=args.field, quiet=args.quiet)
    t = time.time() - t
    if not args.quiet: print("\n All done. Total ellapsed time: %.2f\n" %t)
