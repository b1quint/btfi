#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    imarith.py
    
    This file contains methods and classes used to act exactly as the 'imarith'
    task used on IRAF.
    
    by Bruno C. Quint
    September 2013
"""

from __future__ import division, print_function;

about = "Perform arithmetic operation on 2 input fits files creating " + \
        "a new output file. Operators are '+', '-', '*', '/' and input " + \
        "may be data-cube/data-cube, data-cube/image or image/image. " + \
        "Lists of files are also accepted if the first character is " + \
        "preceded by a '@'.";   

def check_list_of_files_size(_file1_, _file2_):
    """
    Method that checks if two list of files _file1_ and _file2_ have the same 
    number of files inside it.
    
    @param _file1_:    String that contains the path to the _file_ that will be checked.
    @param _file2_:    String that contains the path to the _file_ that will be checked.
    """
    from numpy import loadtxt;
    
    if is_list_of_files(_file1_) and is_list_of_files(_file2_):
        _file1_ = loadtxt(_file1_[1:], dtype=str);
        _file2_ = loadtxt(_file2_[1:], dtype=str);
        assert _file1_.size == _file2_.size;
    
    return;

def load_list(list_name):
    """
    Load a text file containing a list of files. The first char shall be an '@'.
    
    @param list_name:     The name of the list file.
    """
    from numpy import loadtxt;
    
    assert (list_name[0] == '@');
    
    return loadtxt(list_name[1:], dtype=str);

def get_fits_data(fits_name):
    """
    Simple method to handle IO error while reading FITS files.
    
    @param fits_name:     Name of the fits file.
    """
    from pyfits import getdata;
    try:
        foo = getdata(fits_name);
    except IOError:
        from sys import exit;
        print("\n '%s' is not a valid FITS file." % fits_name);
        print(" Please check it and run again.");
        print(" Leaving now.\n");
        exit();
    return foo;
    
def imarith(file1, operator, file2, output, verbose=True, keep=1):
    """
    @param file1:       Input fits image name 1.
    @param operator:    Operator ['+','-','*','/'].
    @param file2:       Input fits image name 2.
    @param output:      Output image.
    @keyword verbose:    Turn verbose on? [True]/False.
    @keyword keep:      Choose whether keep header from [file1] or file2.
    """
    from pyfits import getdata, getheader, writeto;
    from time import strftime;
    
    # Reading data
    data1 = get_fits_data(file1);
    data2 = get_fits_data(file2);
    header = getheader(file2) if keep is 2 else getheader(file1);
    
    # Applying operation
    if verbose: print(" %s %s %s = %s" % (file1, operator, file2, output));
    
    
    if (data1.ndim is 3) and (data2.ndim is 2):
        data2 = data2.reshape((1, data2.shape[0], data2.shape[1]));
        data2 = data2.repeat(data1.shape[0], axis=0);
    if (data1.ndim is 3) and (data2.ndim is 1):
        data2 = data2.reshape((data2.shape[0], 1, 1));
        data2 = data2.repeat(data1.shape[1], axis=1);
        data2 = data2.repeat(data1.shape[2], axis=2);
    
    assert data1.shape == data2.shape; 
    if operator is '+': data1 = data1 + data2;
    elif operator is '-': data1 = data1 - data2;
    elif operator is '*': data1 = data1 * data2;
    elif operator is '/': data1 = data1 / data2;
    else: raise (ValueError, "Invalid operator.");
    
    # Updating header
    try:
        header.set('', '');
    except AttributeError:
        from sys import exit;
        from pyfits import __version__
        print(" Sorry, your PyFits version %s is not supported."% __version__);
        print(" Please, download PyFits v3.0 or superior." );
        print(" Leaving now.\n " );
        exit();
        
    header.set('HISTORY', 'imarith applied on %s' % strftime('%Y-%m-%d %H:%M:%S %Z'));
    header.set('HISTORY', file1, 'First file used on imarith.py.');
    header.set('HISTORY', operator, 'Operator used on imarith.py.');
    header.set('HISTORY', file2, 'Second file used on imarith.py.');
    
    # Writing output
    output = safesave(output);
    writeto(output, data1, header);
    
    return None;

def is_list_of_files(_file_):
    """
    Method that checks if _file_ is, actually, a text file containing
    a set of file names.
    
    @param _file_:    String that contains the path to the _file_ that will be checked.
    """
    from os.path import basename;
    
    _file_ = basename(_file_);
    _file_ = True if _file_[0] == '@' else False;
    
    return _file_;

def safesave(name, overwrite=None, verbose=False):
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
    if v: print(" Writing to output file %s" %name);

    while os.path.exists(name):

        if overwrite in ['y', 'Y', True]:
            if v or verbose: print(" Overwriting %s file." %name);
            os.remove(name);

        elif overwrite in ['', 'n', 'N', False]:
            name = raw_input(" Please, enter a new filename:\n   > ");

        elif overwrite in ['q']:
            if v: print(" Exiting program.");
            sys.exit();

        else:
            overwrite = raw_input("   '%s' file exist. Overwrite? (y/[n])"%name);
            if v: print("   Writing data-cube to %s" %name);

    return name

if __name__ == '__main__':
    
    import argparse;
    import time;
    
    parser = argparse.ArgumentParser(description=about);
    parser.add_argument('input1', metavar='input1', type=str, help="First input FITS file or list of files to be used on the operation.");
    parser.add_argument('operator', metavar='operator', type=str, help="Operator.");
    parser.add_argument('input2', metavar='input2', type=str, help="Second input FITS file list of files to be used on the operation.");
    parser.add_argument('output', metavar='output', type=str, help="Output FITS file or list of files used on the operation.");
    parser.add_argument('-q','--quiet', action='store_true', help="Run quietly.");
    args = parser.parse_args();
    
    v = not args.quiet;
    if v: t = time.time();
    if v: print("\n imarith.py\n Adapted from IRAF's routine with the same name.");
    if v: print(" by Bruno C. Quint - Sep 2013");
    
    if is_list_of_files(args.input1):
        assert is_list_of_files(args.output);
        check_list_of_files_size(args.input1, args.output);
        list1 = load_list(args.input1);
        list3 = load_list(args.output);
        
        if is_list_of_files(args.input2):
            check_list_of_files_size(args.input1, args.input2);
            list2 = load_list(args.input1);
            for i in range(list1.size):
                imarith(list1[i], args.operator, list2[i], list3[i], verbose=v);
        
        else:
            for i in range(list1.size):
                imarith(list1[i], args.operator, args.input2, list3[i], verbose=v);
    else:    
        imarith(args.input1, args.operator, args.input2, args.output, verbose=v);
    
    if v: t = time.time() - t;
    if v: print(" Total elapsed time: %0d:%0d:%0.1f" % (t // 3600, t // 60, t % 60));
    if v: print(" All done.\n");
    
