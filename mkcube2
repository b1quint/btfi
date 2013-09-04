#!/usr/bin/python
#!-*-coding:utf8-*-
"""
    Build Cube v2.0.1

    This script will join the images given as argument to it in a single
    data-cube. The order of the frames within the data-cube will follow the
    alpha-numerical order of the input files.

    usage:
        python % input1.fits input2.fits ... inputN.fits


    v2.0.1 - 'memmap' is now an option.
           - 'os' was not being used. It was removed.

    by Quint, B. C.
    May 2012
"""

import argparse
import numpy
import pyfits
import scipy.signal
import sys

#==============================================================================
# Method Safe Save
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

    v = False if (overwrite is True) else True
    if v: print "  Writing to output file %s" %name

    while os.path.exists(name):

        if overwrite in ['y', 'Y', True]:
            if v or verbose: print "   Overwriting %s file." %name
            os.remove(name)

        elif overwrite in ['', 'n', 'N', False]:
            name = raw_input("   Please, enter a new filename:\n   > ")

        elif overwrite in ['q']:
            if v: print "   Exiting program."
            sys.exit()

        else:
            overwrite = raw_input("   '%s' file exist. Overwrite? (y/[n])"%name)
            if v: print "   Writing data-cube to %s" %name

    return name

#==============================================================================
# Main Thread
if __name__ == '__main__':

    #--------------------------------------------------------------------------
    # Parsing Arguments
    parser = argparse.ArgumentParser(description="Build a data-cube from image files.")
    parser.add_argument('-b','--bin', metavar='bin', type=int, default=1, help="bin size in spacial domains.")
    parser.add_argument('-i','--ipf', metavar='ipf', type=int, default=1, help="number of images per frame.")
    parser.add_argument('-x','--xlim', metavar='xlim', type=int, nargs=2, default=[0, None], help="x limits to be used.")
    parser.add_argument('-y','--ylim', metavar='ylim', type=int, nargs=2, default=[0, None], help="y limits to be used.")
    parser.add_argument('-z','--zlim', metavar='zlim', type=int, nargs=2, default=[0, None], help="z limits to be used.")
    parser.add_argument('-o','--output', metavar='output', type=str, default="cube.fits", help="Name of the output cube.")
    parser.add_argument('-s','--skip', metavar='skip', type=int, default=0, help="# of frames to skip at the beggining of the cube")
    parser.add_argument('-a','--algorithm', metavar='algorithm', type=str, default='average', help="Algorithm used when combining images per frame (average | median)")
    parser.add_argument('-m','--memmap', action='store_true', help="Open/create fits files using MEMMAP.")
    parser.add_argument('-q','--quiet', action='store_true', help="Run quietly.")
    parser.add_argument('files', metavar='files', type=str, nargs='+', help="input filenames.")

    args = parser.parse_args()


    #--------------------------------------------------------------------------
    # Parsing Arguments
    v = not args.quiet
    if v: print "  "
    if v: print "  Build Data-Cube.py"
    if v: print "  by Bruno C. Quint"
    if v: print "  Mar 2013 - Version 2.0.1 "
    if v: print "  "
    if v: print "  Starting program."

    listOfFiles   = args.files
    listOfFiles   = sorted(listOfFiles)
    numberOfFiles = len(listOfFiles)
    if v: print "  %s files will be used." % numberOfFiles
    if v: print "  Ok."
#
#    #--------------------------------------------------------------------------
#    # Setting options
    binStep = args.bin
    imagesPerFrame = args.ipf
    if v: print "  "
    if v: print "  XY binning is: %d" % binStep
    if v: print "  It will be used %d images per frame." % imagesPerFrame

    #--------------------------------------------------------------------------
    # Allocating Memory
    if v: print "  "
    if v: print "  Alocating memory..."
    header   = pyfits.getheader(listOfFiles[0])
    xmin, xmax = args.xlim
    ymin, ymax = args.ylim
    zmin, zmax = args.zlim
    imWidth  = (header.get("NAXIS1") / binStep) if args.xlim==[0, None] else (xmax - xmin) / binStep
    imHeight = (header.get("NAXIS2") / binStep) if args.ylim==[0, None] else (ymax - ymin) / binStep
    imDepth  = numberOfFiles / imagesPerFrame if args.zlim==[0, None] else (zmax - zmin)
    mask  = numpy.ones((binStep,binStep))
    if v: print "  Dimensions of the cube: "
    if v: print "  %d pixels wide" % imWidth
    if v: print "  %d pixels high" % imHeight
    if v: print "  %d frames in z" % imDepth

    if args.algorithm == "average": combFunction = numpy.average
    elif args.algorithm == "median": combFunction = numpy.median
    else: raise IOError, "Wrong algorithm for combining images per frame. Available options are 'average' or 'median'."

    try: cube     = numpy.zeros((imDepth, imHeight, imWidth), dtype=numpy.float32)
    except ValueError, err:
        sys.stderr.write("  Sorry. The cube to be built is too big. Try making a smaller one.\n\n")
        sys.exit()

    if v: print "  Ok."

#    #--------------------------------------------------------------------------
#    # Building Cube
    if v: print "  Building data-cube..."
    if imagesPerFrame == 1:
        for i in range(imDepth):
            current_file = listOfFiles[i + zmin]
            try: frame = pyfits.getdata(current_file, memmap=args.memmap)
            except KeyError: frame = pyfits.getdata(current_file)
            for s in range(binStep):
                cube[i,:,:] += frame[s+ymin:ymax:binStep,s+xmin:xmax:binStep]
            if v:
                sys.stdout.write("\r    %d%%" % (100*(i+1)/imDepth))
                sys.stdout.flush()
    else:
        for i in range(imDepth):
            dummy_cube = numpy.zeros((imagesPerFrame, header.get("NAXIS2") , header.get("NAXIS1")), dtype=float)
            for j in range(imagesPerFrame):
                current_file = listOfFiles[imagesPerFrame * (i + zmin) + j]
                try: frame = pyfits.getdata(current_file, memmap=args.memmap)
                except KeyError: frame = pyfits.getdata(current_file)
                dummy_cube[j,:,:] = frame
            frame = combFunction(dummy_cube, axis=0)
            frame = frame[ymin:ymax,xmin:xmax]
            del dummy_cube

            if binStep != 1:
                frame = scipy.signal.convolve2d(frame, mask,'same')[1::binStep,1::binStep]
            cube[i,:,:] = frame
            del frame

            if v:
                sys.stdout.write("\r    %d%%" % (100*(i+1)/imDepth))
                sys.stdout.flush()
    if v: print "  Ok.\n"

    #--------------------------------------------------------------------------
    # Fixing header
    header.update('CRVAL1', xmin+1)
    header.update('CRVAL2', ymin+1)
    header.update('CRVAL3', zmin+1)

    header.update('CDELT1', binStep)
    header.update('CDELT2', binStep)
    header.update('CDELT3', 1)


    #--------------------------------------------------------------------------
    # Writing file
    cubename = safesave(args.output)
    pyfits.writeto(cubename, cube, header)
    del cube
    del header
    if v: print "  All done.\n"


