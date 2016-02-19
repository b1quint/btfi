#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
     Phase-Map Apply - A script to apply a phase-map on a data-cube.
     by Bruno Quint (bquint@astro.iag.usp.br)
     and Fabricio Ferrari (fabricio@ferrari.pro.br)
     version 0.0 - Feb 2014
"""

from __future__ import division, print_function

import argparse
import astropy.io.fits as pyfits
import numpy
import os
import sys
import time

from scipy.interpolate import UnivariateSpline

def main():

    # Setting Options ---------------------------------------------------------
    parser = argparse.ArgumentParser(description="Apply a phase-map on" + \
                                     "a data-cube.")

    parser.add_argument('-o', '--output', metavar='output', type=str,
                        default=None, help="Name of the output corrected cube")

    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run it quietly.")

    parser.add_argument('-n', '--npoints', type=int, default=10,
                        help="Number of points in the re-sampling for channel [10].")

    parser.add_argument('cube_file', metavar='cube_file', type=str,
                        help="Input calibration cube filename.")

    parser.add_argument('map_file', metavar='map_file', type=str,
                        help="Input phase-map image filename.")

    args = parser.parse_args()
    v = not args.quiet
    loading = [' ','-','\\','|','/']

    # Printing program header --------------------------------------------------
    if v:
        start = time.time()
        print("\n Phase-Map Apply")
        print(" by Bruno Quint & Fabricio Ferrari")
        print(" version 0.0 - Feb 2014")

    root_dir = os.path.dirname(args.cube_file)
    cube_file = args.cube_file
    map_file = args.map_file

    if args.output is None:
        out_file = 'phc_' + os.path.split(args.cube_file)[-1]
    else:
        out_file = args.output

    if v:
        print(" \n Root dir: %s" % root_dir)
        print(" Cube to be corrected: %s" % cube_file)
        print(" Phase-map to be applied: %s" % map_file)
        print(" Output corrected cube: %s" % out_file)

    # Reading input data ------------------------------------------------------
    if v:
        print("\n Reading cube to be corrected.")

    data_cube = pyfits.open(cube_file)[0]

    if v:
        print(" Done.")
        print("\n Reading phase-map to be applied.")

    phase_map = pyfits.open(map_file)[0]

    if v:
        print(" Done.")

    # Checking data -----------------------------------------------------------
    if data_cube.data[0].shape != phase_map.shape:
        print("[!] Cube and map does not have matching width and height.")
        print("[!] Leaving now.\n")
        sys.exit()

    if data_cube.data.ndim != 3:
        print("[!] Cube file is not really a cube.")
        print("[!] Leaving now.\n")
        sys.exit()

    if phase_map.data.ndim != 2:
        print("[!] Map file is not really an image.")
        print("[!] Leaving now.\n")
        sys.exit()

    check_instrument(cube_file)
    mode = check_mode(cube_file)

    ## Phase-Correction for iBTF data-cube ------------------------------------
    if mode.lower() in ['ibtf']:

        # Padding data-cube with zeros ----------------------------------------
        if v:
            print(" Padding data-cube with zeros")

        phm_max = round(abs(phase_map.data).max()) + 1
        L, M, N = data_cube.data.shape
        pad = numpy.zeros((phm_max, M, N))

        if v:
            print(" Cube shape before paddding: %d x %d x %d" % (N, M, L))
            print(" %d frames will be added." % (2 * phm_max))

        data_cube.data = numpy.vstack((pad, data_cube.data))
        data_cube.data = numpy.vstack((data_cube.data, pad))
        L, M, N = data_cube.data.shape

        try:
            data_cube.header['CRPIX3'] = data_cube.header['CRPIX3'] + phm_max
        except KeyError:
            data_cube.header['CRPIX3'] = L

        if v:
            print(" Cube shape after padding: %d x %d x %d" % (N, M, L));
            print(" Done.");

        # Applying phase-map --------------------------------------------------
        if v:
            print("\n Applying phasemap")
        for i in range(M):
            for j in range(N):
                if v:
                    temp = (((i + 1) * 100.00 / M))
                    sys.stdout.write('\r  %2d%% ' % temp)
                    sys.stdout.write(loading[int(temp * 10 % 5)])
                    sys.stdout.flush()
                spec = data_cube.data[:,i,j]
                shift = phase_map.data[i,j]
                data_cube.data[:,i,j] = shiftSpectrum(spec, shift, args.npoints)
        if v: print(" Done.")

    ## Phase-Correction for Fabry-Perot data-cube ------------------------------
    elif mode.lower() in ['fabry-perot', 'fp']:

        M = data_cube.header['NAXIS1']
        N = data_cube.header['NAXIS2']

        ref_x = phase_map.header['PHMREFX']
        ref_y = phase_map.header['PHMREFY']
        units = phase_map.header['PHMUNIT']
        sample = float(phase_map.header['PHMSAMP'])

        # Reading the Free-Spectral-Range --------------------------------------
        try:
            if v:
                print(" Reading free-spectral-range from cube header.")
            # TODO add an option to use the FSR found while extracting
            # TODO the phase-map or while fitting it.
            # TODO or even to give the option for the user to enter it.
            # FSR = phase_map.header['PHMFITSR']
            FSR = phase_map.header['PHMFSR']
            if v:
                print(" Free Spectral Range = %.2f %s" % (FSR, units))

        except (KeyError):
            print(" Please, enter the free-spectral-range in %s units" % units)
            FSR = input(" > ")

        FSR = round(FSR / sample) # From BCV to Channels
        if v:
            print(" Free-Spectral-Range is %d channels" % FSR)

        fsr = FSR * args.npoints # From Channels to nPoints
        fsr = int(round(fsr))
        if v:
            print(" Free-Spectral-Range is %d points" % fsr)

        # Assure that the reference spectrum will not be moved ----------------
        try:
            phase_map.data = phase_map.data - phase_map.data[ref_y, ref_x]
        except IndexError:
            print("[!] Reference pixel out of field.")
            print("[!] Skipping reference pixel map subtraction.")
            pass
        phase_map.data = -1 * phase_map.data

        # Converting phase-map values to channels ------------------------------
        phase_map.data = phase_map.data / sample

        # Converting phase-map from channels to number of points --------------
        phase_map.data = phase_map.data * args.npoints

        # Applying phase-map --------------------------------------------------
        if v:
            print("\n Applying phase-map:")

        z = numpy.arange(data_cube.header['NAXIS3'])
        new_z = numpy.arange(0, z.size, 1.0 / args.npoints)

        for i in range(M):
            for j in range(N):

                # Extracting a spectrum
                spec = data_cube.data[:,j,i]
                dz = phase_map.data[j,i]

                # Re-sample spectrum
                spline = UnivariateSpline(z, spec, s=0.0)
                new_spec = spline(new_z)

                # Cutting spectrum
                new_z = new_z[0:fsr+1]
                new_spec = new_spec[0:fsr+1]

                # Shifting spectrum
                new_spec = numpy.roll(new_spec, int(dz))

                # Under-sampling spectrum
                spline = UnivariateSpline(new_z, new_spec, s=0.0)
                spec = spline(z)

                # Storing new spectrum
                data_cube.data[:,j,i] = spec

                # Giving a feedback to the user
                if v:
                    temp = (((i + 1) * 100.00 / M))
                    sys.stdout.write('\r  %2.2f%% ' % temp)
                    sys.stdout.flush()

        end_of_cube = min(int(round(FSR)), data_cube.data.shape[0])
        data_cube.data = data_cube.data[0:end_of_cube, :, :]
        if v: print(" Done.")

    else:
        sys.exit()

    # Saving corrected data-cube ----------------------------------------------
    if v:
        print("\n Writing output to file %s." % out_file);
    data_cube.writeto(out_file, clobber=True)
    if v:
        print(" Done.");
        end = time.time() - start
        print("\n Total time ellapsed: %02d:%02d:%02d" % \
              (end // 3600, end % 3600 // 60, end % 60));
        print(" All done!\n");

def check_instrument(filename, instrument='btfi', keyword='INSTRUME'):
    """
    Method written to check the instrument.
    """
    header = pyfits.getheader(filename)

    # First check if the keyword exists
    if not keyword in header:
        print("")
        print(" Instrument type not recognized.")
        print(" Do you want to proceed? [Y or n]")

        answer = '.'
        while answer.lower() not in ' yn':
            answer = raw_input('? ')

        if answer.lower() == 'n':
            print(" Leaving now.\n")
            sys.exit()
        else:
            return

def check_mode(filename, keyword='INSTRMOD'):
    """
    Return if BTFI was obtained with a Fabry-Perot or with the iBTF.
    """
    from astropy.io import fits as pyfits

    header = pyfits.getheader(filename)

    if keyword not in header:
        warning("Instrument mode not found.")
        instrument_mode = ''
        while instrument_mode.lower() not in ['ibtf', 'fp']:
            instrument_mode = raw_input(" Enter 'ibtf' or 'fp': ")
    else:
        if header[keyword].upper() in ['IBTF']:
            instrument_mode = 'ibtf'

        if header[keyword].upper() in ['FP', 'FABRY-PEROT']:
            instrument_mode = 'fp'

    return instrument_mode

## Method shiftSpectrum ========================================================
def shiftSpectrum(spec, dz, nPoints=100):

    dzSign = -numpy.sign(dz)
    dz = abs(dz)
    dzPoints = int(dz * nPoints)
    if dzPoints is 0: return spec

    # Get the spectrum from cube
    z = numpy.arange(spec.size)
    spline = UnivariateSpline(z, spec, s=0.0)

    # Re-sample spectrum
    newZ    = numpy.linspace(z[0], z[-1], z.size * nPoints)
    newSpec = spline(newZ)

    # Add padded borders
    newSpec = numpy.append(numpy.zeros(dzPoints), newSpec)
    newSpec = numpy.append(newSpec, numpy.zeros(dzPoints))

    # Shifting spectrum
    newSpec = numpy.roll(newSpec, int(dzSign * dzPoints))

    # Cutting Spectrum
    newSpec = newSpec[dzPoints:-dzPoints]

    # Down-sampling
    spline = UnivariateSpline(newZ, newSpec, s=0.0)
    spec   = spline(z)

    return spec

## Method shift_spectrum ========================================================
def shift_spectrum(spec, dz, fsr=-1, sample=1.0, n_points=100):
    """
    @param spec: 1D numpy.array representing a spectrum to be shifted.
    @param dz: how big is the shifting.
    @keyword fsr: a float representing the free-spectra-range in sample units.
    @keyword sample: a float representing the increment between each channel.
    @keyword n_points: number of points that will be used for super-sampling.
    """

    dzSign = -numpy.sign(dz)
    dz = abs(dz) / sample # From cube units to channels
    dzPoints = int(dz * n_points) # From channels to new sample units

    index = fsr / sample # From cube units to channels
    index = index * n_points # From channels to new sample units

    if dzPoints is 0:
        return spec

    # Get the spectrum from cube
    z = numpy.arange(spec.size)
    spline = UnivariateSpline(z, spec, s=0.0)

    # Re-sample spectrum
    newZ = numpy.linspace(z[0], z[-1], z.size * n_points)
    newSpec = spline(newZ)

    # Cutting Spectrum
    newSpec = newSpec[0:fsr]

    # Shifting spectrum
    newSpec = numpy.roll(newSpec, int(dzSign * dzPoints))

    # Down-sampling
    spline = UnivariateSpline(newZ, newSpec, s=0.0)
    spec   = spline(z)

    return spec

def error(my_string):
    s = bcolors.FAIL + '[ERROR] ' + bcolors.ENDC
    s = s + my_string
    print(s)
    return

def warning(my_string):
    s = bcolors.WARNING + '[WARNING] ' + bcolors.ENDC
    s = s + my_string
    print(s)
    return

class bcolors:

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

#===============================================================================
if __name__ == '__main__':
    main()
