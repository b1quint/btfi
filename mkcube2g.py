#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    Build Cube v2g

    This script will join the images given as argument to it in a single
    data-cube. The order of the frames within the data-cube will follow the
    alpha-numerical order of the input files.

    v2g - Changed header.update to header['xxx'] following astropy.io.fits
          directives.
    v2f - Error handling when finding a corrupted FITS file.
          (KeyError, IOError, IndexError)
    v2e - Print the name of the algorithm used for combining the images.
        - Added an option to subtract BIAS while combining.
        - Added an option to subtract DARK while combining.
        - Added an option to divide by FLAT while combining.
        - Identify whether the data-cube is a FP data-cube or an iBTF data-cube.
        - Added exception for corrupted files. Filled with zeros.
    v2d - A 'main' function was created to organize better the code.
        - Code reorganized according to PEP 8.
        - Inserted 'from __future__' statement.
    v2c - 'memmap' is now an option.
        - 'os' was not being used. It was removed.

    by Quint, B. C.
    Jan 2014

    TODO:
    - Organize the data-cube using cards in the header.
    - If 2D Flat, divide all frames by it. If 3D Flat, match each frame.
"""

from __future__ import division, print_function

def main():
    """
    This is actually what happens when the program run.
    """

    import argparse
    import astropy.io.fits as pyfits
    import numpy
    import scipy.signal
    import sys

    # Parsing Arguments -------------------------------------------------------
    parser = argparse.ArgumentParser(description="Build a data-cube from image files.")
    parser.add_argument('-a','--algorithm', metavar='algorithm', type=str, default='average', help="Algorithm used when combining images per frame (average | median | sum)")
    parser.add_argument('-b','--bin', metavar='bin', type=int, default=1, help="bin size in spacial domains.")
    parser.add_argument('-B','--bias', default=None, type=str, help='Bias image to be subtracted on each input image.')
    parser.add_argument('-d','--dark', default=None, type=str, help='Dark image to be subtracted on each input image.')
    parser.add_argument('-f','--flat', default=None, type=str, help='Flat image for normalization')
    parser.add_argument('-i','--ipc', metavar='ipc', type=int, default=1, help="number of images per channel.")
    parser.add_argument('-m','--memmap', action='store_true', help="Open/create fits files using MEMMAP.")
    parser.add_argument('-o','--output', metavar='output', type=str, default="cube.fits", help="Name of the output cube.")
    parser.add_argument('-q','--quiet', action='store_true', help="Run quietly.")
    parser.add_argument('-s','--skip', metavar='skip', type=int, default=0, help="# of frames to skip at the beggining of the cube")
    parser.add_argument('-x','--xlim', metavar='xlim', type=int, nargs=2, default=[0, None], help="x limits to be used.")
    parser.add_argument('-y','--ylim', metavar='ylim', type=int, nargs=2, default=[0, None], help="y limits to be used.")
    parser.add_argument('-z','--zlim', metavar='zlim', type=int, nargs=2, default=[0, None], help="z limits to be used.")
    parser.add_argument('files', metavar='files', type=str, nargs='+', help="input filenames.")

    args = parser.parse_args()

    v = not args.quiet
    if v:
        print("\n mkCube")
        print(" by Bruno Quint (bquint@astro.iag.usp.br)")
        print(" Apr 2014 - Version 2g")
        print("\n Starting program.")

    list_of_files = args.files
    list_of_files = sorted(list_of_files)
    number_of_files = len(list_of_files)

    header = pyfits.getheader(list_of_files[0])
    if v:
        try:
            print(" Instrument mode: %s" % header['INSTRMOD'].upper())
        except KeyError:
            print(" Instrument mode not found. No calibration for third axis will be added.")
        print(" %s files will be used." % number_of_files)
        print(" Ok.")

    bin_step = args.bin
    images_per_channel = args.ipc
    if v:
        print("\n XY binning is: %d" % bin_step)
        print(" It will be used %d images per channel." % images_per_channel)

    # Allocating Memory -------------------------------------------------------
    if v:
        print("\n Allocating memory.")

    xmin, xmax = args.xlim
    ymin, ymax = args.ylim
    zmin, zmax = args.zlim

    if args.xlim == [0, None]:
        imWidth = (header.get("NAXIS1") // bin_step)
    else:
        imWidth = (xmax - xmin) // bin_step

    if args.ylim == [0, None]:
        imHeight = (header.get("NAXIS2") // bin_step)
    else:
        imHeight = (ymax - ymin) // bin_step

    if args.zlim == [0, None]:
        imDepth = number_of_files // images_per_channel
    else:
        imDepth = (zmax - zmin)

    mask = numpy.ones((bin_step,bin_step))
    if v:
        print(" Dimensions of the cube: ")
        print(" %d pixels wide" % imWidth)
        print(" %d pixels high" % imHeight)
        print(" %d frames in z" % imDepth)
        print("\n Combining images on each frame using algorithm: %s" %
              args.algorithm)

    # print how much RAM we will use
    memory_needed_mb = imWidth * imHeight * imDepth * 4 / ( 1024 * 1024 )
    print ( " We need %d MB of RAM.\n" % memory_needed_mb )

    if args.algorithm == "average":
        combFunction = numpy.average
    elif args.algorithm == "median":
        combFunction = numpy.median
    elif args.algorithm == "sum":
        combFunction = numpy.sum
    else:
        raise IOError, "Wrong algorithm for combining images" + \
            "per frame. Available options are 'average' or 'median'."

    try:
        cube = numpy.zeros((imDepth, imHeight, imWidth), dtype=numpy.float32)
    except ValueError:
        sys.stderr.write(" Sorry. The cube to be built is too big.")
        sys.stderr.write(" Try making a smaller one.\n\n")
        sys.exit()

    except MemoryError:
        sys.stderr.write ( " ERROR: Numpy threw a MemoryError exception.\n" )
        sys.stderr.write ( "        It is probable that the OS could not allocate a large enough chunk of memory.\n\n" )
        sys.exit ( )

    if v:
        print(" Ok.")

    # Loading bias file -------------------------------------------------------
    if args.bias != None:
        print(' Loading "%s" file for BIAS subtraction.' % args.bias)
        bias = pyfits.getdata(args.bias)
        header.set('BIASFILE', args.bias, 'File used for BIAS subtraction.')
    else:
        bias = numpy.zeros((imHeight, imWidth))

    # Loading dark file -------------------------------------------------------
    if args.dark != None:
        print(' Loading "%s" file for DARK subtraction.' % args.dark)
        dark = pyfits.getdata(args.dark)
        header.set('DARKFILE', args.bias, 'File used for DARK subtraction.')
    else:
        dark = numpy.zeros((imHeight, imWidth))

    # Loading flat file -------------------------------------------------------
    if args.flat != None:
        print(' Loading "%s" file for FLAT normalization.' % args.flat)
        flat = pyfits.getdata(args.flat)
        header.set('FLATFILE', args.flat, 'File used for FLAT normalization.')
    else:
        flat = numpy.ones((imHeight, imWidth))

    # Building Cube ------------------------------------------------------------

    if v:
        print("\n Building data-cube.")

    if images_per_channel == 1:
        for i in range(imDepth):
            current_file = list_of_files[i + zmin]
            try:
                frame = pyfits.getdata(current_file, memmap=args.memmap)
            except (KeyError, IOError, IndexError):
                print(" Warning: %s file corrupted." % current_file)
                frame = numpy.zeros((header['NAXIS2'], header['NAXIS1']))
            for s in range(bin_step):
                frame = frame - bias
                frame = frame - dark
                frame = frame / flat
                cube[i,:,:] += frame[s+ymin:ymax:bin_step,s+xmin:xmax:bin_step]
            sys.stdout.write("\r    %d%%" % (100*(i+1)/imDepth))
            sys.stdout.flush()
    else:
        for i in range(imDepth):
            dummy_cube = numpy.zeros((1,
                                      header.get("NAXIS2"),
                                      header.get("NAXIS1")),
                                      dtype=float)
            for j in range(images_per_channel):
                current_file = list_of_files[images_per_channel * (i + zmin) + j]
                try:
                    frame = pyfits.getdata(current_file, memmap=args.memmap)
                    frame = frame - bias
                    frame = frame - dark
                    frame = frame / flat
                    frame = numpy.reshape(frame, (1, frame.shape[0], frame.shape[1]))
                    dummy_cube = numpy.append(dummy_cube, frame, axis=0)
                except (KeyError, IOError, IndexError):
                    print(" Warning: %s file corrupted." % current_file)
            frame = combFunction(dummy_cube, axis=0)
            frame = frame[ymin:ymax,xmin:xmax]
            del dummy_cube

            if bin_step != 1:
                frame = scipy.signal.convolve2d(frame, mask,'same')[1::bin_step,1::bin_step]
            cube[i,:,:] = frame
            del frame

            if v:
                sys.stdout.write("\r    %d%%" % (100*(i+1)/imDepth))
                sys.stdout.flush()
    if v:
        print(" Ok.\n")

    # Fixing header -----------------------------------------------------------
    header['CRPIX1'] = 1
    header['CRVAL1'] = xmin + 1
    header['CDELT1'] = bin_step
    header['CTYPE1'] = 'LINEAR'
    header['CUNIT1'] = 'PIXEL'

    header['CRPIX2'] = 1
    header['CRVAL2'] = ymin + 1
    header['CDELT2'] = bin_step
    header['CTYPE2'] = 'LINEAR'
    header['CUNIT2'] = 'PIXEL'

    if header['INSTRMOD'].lower() in ['fp', 'fabry-perot']:
        header['CRPIX3'] = 1
        header['CRVAL3'] = float(header['FPZINIT'])
        header['CDELT3'] = float(header['FPZDELT'])
        header['CTYPE3'] = 'LINEAR'
        header['CUNIT3'] = 'BNV'
        header['C3_3'] = float(header['FPZDELT'])

    elif header['INSTRMOD'].lower() in ['ibtf']:
        header['CRPIX3'] = 1
        header['CRVAL3'] = float(header['TFAINIT'])
        header['CDELT3'] = float(header['TFADELT'])
        header['CTYPE3'] = 'LINEAR'
        header['CUNIT3'] = 'degrees'
        header['C3_3'] = float(header['TFADELT'])

    else:
        print("[!] Invalid BTFI Instrument Mode.")
        print("[!] Dummy calibration will be added to the data-cube.")
        header['CRPIX3'] = 1
        header['CRVAL3'] = 1
        header['CDELT3'] = 1
        header['CTYPE3'] = 'LINEAR'
        header['CUNIT3'] = 'channel'
        header['C3_3'] = 1

    header.set('', '-- mkcube calibration --', before='CRPIX1')

    # Writing file ------------------------------------------------------------
    cubename = safesave(args.output)
    pyfits.writeto(cubename, cube, header)
    del cube
    del header
    if v:
        print("\n All done.\n")

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
    import os
    import sys

    v = False if (overwrite is True) else True
    if v: print("\n Writing to output file %s" % name)

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
            overwrite = raw_input("   '%s' file exist. Overwrite? (y/[n])"%name)
            if v:
                print(" Writing data-cube to %s" %name)

    return name


# Main Thread -----------------------------------------------------------------
if __name__ == '__main__':
    main()

