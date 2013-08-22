#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    hlist.py
    
    This file contains a script used to list all the input FITS files that have
    a FIELD with a value that matches the input.
    
    by Bruno C. Quint (bquint@astro.iag.usp.br)
    August 2013
"""

from __future__ import division, print_function

def hselect(list_of_input_files, header_field, field_value, verbose=False):
    """
    Returns a list with the name of the files contained at 'list_of_files'
    whose 'header_field' as a value matching the 'field_value'.
    
    @param list_of_input_files: A list containing the name of the files that will be analised.
    @param header_field: The name of the FIELD that will be used to select the files.
    @param field_value: The value of the FIELD that will be used to select the files.
    """
    from astropy.io.fits import getheader
    
    import warnings
    warnings.filterwarnings('ignore')
 
    if field_value.isdigit():
        field_value = float(field_value);
 
    selected_files = []
    for _file_ in list_of_input_files:
        try:
            header = getheader(_file_);
            value = header.get(header_field);
            value = value if not value.isdigit() else float(value);  
        except IOError:
            continue
        
        if value == field_value:
            selected_files.append(_file_);
            if verbose: print(_file_);

    return selected_files

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser(description="List all the input FITS files whose FIELD has a value that matches VALUE.");
    parser.add_argument('field', type=str, help="The FIELD used to select the files.");
    parser.add_argument('value', help="The VALUE of the FIELD that serves as criterium to list the files.");
    parser.add_argument('files', type=str, nargs='+', help="Input files to be selected.");
    args = parser.parse_args();
        
    list_of_files = hselect(args.files, args.field, args.value, verbose=True);
