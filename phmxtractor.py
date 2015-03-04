#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    Phase-map Xtractor
    by Bruno C Quint

    v1a - Phase extraction for Fabry-Perot.
    2014.04.16 15:45 - Created an exception for errors while trying to access
                       'CRPIX%' cards on cube's header.
"""
from __future__ import division, print_function

import argparse
import astropy.io.fits as pyfits
import matplotlib.pyplot as pyplot
import numpy
import time
import scipy
import scipy.interpolate as interpolate
import scipy.ndimage as ndimage
import sys

def main():

    # Parse arguments ---------------------------------------------------------
    parser = argparse.ArgumentParser(description="Extracts the phase-map" +
                                     "from a fits file containing a data" +
                                     "-cube.")

    parser.add_argument('-c', '--correlation', action='store_true',
                        help="Use correlation cube? true/[FALSE]")

    parser.add_argument('filename', type=str, help="Input data-cube name.")

    parser.add_argument('-o', '--output', type=str, default=None,
                        help="Name of the output phase-map file.")

    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run program quietly. true/[FALSE]")

    parser.add_argument('-s', '--show', action='store_true',
                        help="Show plots used in the process. true/[FALSE]")

    args = parser.parse_args()

    # Starting program --------------------------------------------------------
    v = not args.quiet
    start = time.time()

    if v:
        print("")
        print(" Phase-Map Extractor")
        print(" by Bruno Quint & Fabricio Ferrari")
        print(" version 0.1c - May 2014")
        print(" Extracting phase-map from file: %s" % args.filename)

    # Checking input data -----------------------------------------------------
    if v:
        print(" Checking data-cube for phase-correction.")
    check_dimensions(args.filename)
    check_instrument(args.filename)

    # Selecting BTFI mode and extracting phase-map -----------------------------
    mode = check_mode(args.filename)
    if mode == 'ibtf':
        PhaseMap_iBTF(args.filename, correlation=args.correlation,
                      show=args.show, verbose=v)
    elif mode == 'fp':
        PhaseMap_FP(args.filename, correlation=args.correlation,
                    show=args.show, verbose=v)

    # All done! ---------------------------------------------------------------
    end = time.time() - start
    if v:
        print("\n  Total time ellapsed: %02d:%02d:%02d" %
              (end // 3600, end % 3600 // 60, end % 60))
        print("  All done!\n")

def check_dimensions(filename, dimensions=3, keyword='NAXIS'):
    """
    Method written to check the dimensions of the input fits data.
    """
    header = pyfits.getheader(filename)

    if not 'NAXIS' in header:
        data = pyfits.getdata(args.filename)
        ndim = data.ndim
    else:
        ndim = header['NAXIS']

    if ndim is not dimensions:
        print(" INVALID DATA: wrong number of dimensions")
        print(" Leaving now!")
        print("")
        sys.exit()
    else:
        return

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

    # Then check if it is the right instrument
    if header[keyword].lower() is not instrument:
        print(" Wrong instrument. Proceed? [Y or n]")

        answer = '.'
        while answer.lower() not in ' yn':
            answer = raw_input('? ')

        if answer.lower() == 'n':
            print(" Leaving now.\n")
            sys.exit()

    return

def check_mode(filename, keyword='INSTRMOD'):
    """
    Return if BTFI was obtained with a Fabry-Perot or with the iBTF.
    """
    header = pyfits.getheader(filename)

    if keyword not in header:
        print("\n Instrument mode not found.")
        instrument_mode = ''
        while instrument_mode.lower() not in ['ibtf', 'fp']:
            instrument_mode = raw_input("? Enter 'ibtf' or 'fp': ")
    else:
        if header[keyword].upper() in ['IBTF']:
            instrument_mode = 'ibtf'

        if header[keyword].upper() in ['FP', 'FABRY-PEROT']:
            instrument_mode = 'fp'

    return instrument_mode

def safe_save(name, extension=None, overwrite=False, verbose=False):
    """
    This is a generic method used to check if a file called 'name' already 
    exists. If so, it starts some interaction with the user.
    
    @param name: the name of the file that will be written in the future.
    
    @keyword extension: check if the file has the following extension. If not,
    it fills for the user. Defaults is None. An example would be 
    extension='.fits'.
    
    @keyword overwrite: if False, this method will interact with the user to 
    ask if 'name' file shall be overwritten or if a new name will be given. If
    True, 'name' file is automatically overwritten.
    
    @keyword verbose: force verbose mode on even when overwrite is automatic. 
    
    v1.0.2 - added 'extension' keyword.
    v1.0.1 - added 'overwrite' keyword.
           - added 'verbose' keyword.
    """
    import os

    if os.path.splitext(name)[1] != extension and extension is not None:
        name = name + extension

    v = False if (overwrite is True) else True
    if v: print('\n Writing to output file "%s"' % name)

    while os.path.exists(name):

        if overwrite in ['y', 'Y', True]:
            if v or verbose:
                print(" Overwriting %s file." % name)
            os.remove(name)

        elif overwrite in ['', 'n', 'N', False]:
            name = raw_input(" Please, enter a new filename:\n > ")
            if os.path.splitext(name)[1] != extension and extension is not None:
                name = name + extension

        elif overwrite in ['q']:
            if v:
                print(" Exiting program.")
            sys.exit()

        else:
            overwrite = raw_input(" '%s' file exist. Overwrite? (y/[n])"%name)
            if v:
                print(" Writing data-cube to %s" %name)

    return name



#==============================================================================
class PhaseMap:

    def __init__(self, filename, **config):

        # Setting main configuration ------------------------------------------
        self.input_file = filename
        self.config = config
        self.verbose = config['verbose']

        if 'show' in config:
            self.show = config['show']
        else:
            self.show = False

        self.loading = [' ','-','\\','|','/']

        # Reading raw data ----------------------------------------------------
        self.print("  Loading data.")
        self.data = pyfits.getdata(filename)
        self.header = pyfits.getheader(filename)
        self.print("  Done.")

        # Reading data-cube configuration -------------------------------------
        self.width = self.header['NAXIS1']
        self.height = self.header['NAXIS2']
        self.depth = self.header['NAXIS3']

        # Reading Z calibration for plotting ----------------------------------
        self.z = self.get_calibration()

        try:
            self.units = self.header['CUNIT3']
        except KeyError:
            self.units = "channels"

        try:
            self.sample = self.header['C3_3']
        except KeyError:
            self.sample = 1.0

        return

    def extract_phase_map(self):
        """
        Extract the phase-map.
        """
        from astropy.io.fits import getdata
        from numpy import argmax, inf, where

        self.print("\n  Starting phase-map extraction.")
        self.print("  Reading data from %s file" % self.extract_from)
        data = getdata(self.extract_from)
        #data = where(data > data.mean() + data.std(), data, -inf)

        phase_map = argmax(data, axis=0) * self.sample
        return phase_map

    def find_reference_pixel(self):
        """Read the reference pixel from header or find it."""
        if ('PHMREFX' in self.header) and ('PHMREFY' in self.header):
            self.print(" \n  Found reference pixel in header.")
            ref_x = self.header['PHMREFX']
            ref_y = self.header['PHMREFY']
            self.print("  Using [%d, %d]" % (self.ref_x, self.ref_y))
        else:
            self.print(" \n  Reference pixel NOT in header.")

            answer = '.'
            while answer.lower() not in ' yn':
                answer = raw_input("? Use central pixel? [Y, n]\n? ")

            if answer.lower() in ' y':
                ref_x = self.width // 2
                ref_y = self.height // 2
            else:
                ref_x = int(raw_input("? Please, enter reference X: "))
                ref_y = int(raw_input("? Please, enter reference Y: "))

            self.print("  Using [%d, %d]" % (ref_x, ref_y))

        return ref_x, ref_y

    def get_calibration(self):
        """
        Return an array with the current calibration.
        """
        z = numpy.arange(self.depth)
        try:
            # The "+ 1" change from fortran like to c like indexing
            z = z - self.header['CRPIX3'] + 1
            z = z * self.header['C3_3']
            z = z + self.header['CRVAL3']

        except KeyError:
            print("! Calibration in third axis not found.")
            print("! I will ignore this step.")

        return z

    def get_fwhm(self):
        """
        Returns the full-width-at-half-maximum.
        """

        # TODO Add an option to choose wheter to use gauss FWHM or measured FWHM

        from scipy.optimize import leastsq
        from scipy.stats import mode

        fsr = self.free_spectral_range / self.header['C3_3']
        z = self.z[:fsr]
        s = self.ref_s[:fsr]
        s = s - mode(s)[0]


        zz = numpy.linspace(z[0], z[-1], 1000)
        ss = interpolate.interp1d(z, s, kind='cubic')
        sss = ss(zz) - ss(zz).max() / 2

        fit_func = lambda p, x: p[0] * numpy.exp(-(x - p[1]) ** 2 / (2 * p[2] ** 2))
        err_func = lambda p, x, y: y - fit_func(p, x)
        p = [sss.max(), zz[sss.argmax()], 10]
        p, _ = leastsq(err_func, p, args=(zz, sss))
        fwhm_gauss = 2.35482 * p[2]

        zzz = zz[sss > 0]
        fwhm_measured = zzz.ptp()

        if self.show:
            pyplot.figure()
            pyplot.title("Measure the FWHM")
            pyplot.plot(z, s, 'bo')
            pyplot.plot(zz, ss(zz), 'b-', lw=2)
            pyplot.plot(zz, sss, 'r-', lw=2, alpha=0.3)
            pyplot.plot(zz, fit_func(p, zz), 'g-', lw=2, alpha=0.3)
            pyplot.axvline(p[1] - fwhm_gauss / 2, ls='--', c='green', lw=2)
            pyplot.axvline(p[1] + fwhm_gauss/ 2, ls='--', c='green', lw=2,
                label='Gauss Fit = %.1f %s' % (fwhm_gauss, self.units))
            pyplot.axvline(p[1] + fwhm_measured/ 2, ls='--', c='red', lw=2)
            pyplot.axvline(p[1] - fwhm_measured/ 2, ls='--', c='red', lw=2,
                label='Definition = %.1f %s' % (fwhm_measured, self.units))
            pyplot.legend(loc='best')
            pyplot.grid()
            pyplot.tight_layout()
            pyplot.show()

        if self.verbose:
            print("")
            print(" Measured FWHM = %.2f %s" % (fwhm_measured, self.units))
            print(" Gauss-fit FWHM = %.2f %s " % (fwhm_gauss, self.units))
            print(" Using the measured FWHM for further calculations.")

        return fwhm_measured

    def get_reference_spectrum(self):
        """
        Get the reference spectrum.
        """
        from scipy.stats import mode

        ref_s = pyfits.getdata(self.input_file)[:,self.ref_y, self.ref_x]
        ref_s = ref_s / ref_s.max()  # Normalize
        ref_s = ref_s - ref_s.mean() # Remove mean to avoid triangular shape
        ref_s = ref_s - mode(ref_s)[0] # Try to put zero on zero

        if self.show:
            pyplot.figure()
            pyplot.title("Reference Spectrum")
            pyplot.plot(self.z, ref_s, 'ko-', label="Reference spectrum")
            pyplot.grid()
            pyplot.xlabel("z [%s]" % self.units)
            pyplot.tight_layout()
            pyplot.show()

        return ref_s

    def get_refx_pixel(self):
        """
        Return the position of the reference X in pixels.
        """

        return

    def print(self, string):
        """
        Print only in verbose mode.
        """
        if self.verbose: print(string)
        return

    def use_correlation(self):
        """
        Use correlation data-cube.
        """
        import numpy

        from astropy.io.fits import getdata, getheader, writeto
        from glob import glob
        from os.path import splitext
        from sys import stdout

        self.print("\n A correlation cube will be used.")
        self.print(" Looking for an existing correlation data-cube in the current folder.")
        candidates = glob("*.fits")

        corr_cube = None
        for candidate in candidates:
            if 'CORRFROM' in getheader(candidate):
                if getheader(candidate)['CORRFROM'] == self.input_file:
                    self.print(" Correlation cube to be used: %s" % candidate)
                    return candidate

        if corr_cube == None:
            self.print(" Correlation cube not found. Creating a new one.")
            data = getdata(self.input_file)
            corr_cube = numpy.empty_like(data)

            x = numpy.arange(self.width)
            y = numpy.arange(self.height)
            X, Y = numpy.meshgrid(x, y)
            x, y = numpy.ravel(X), numpy.ravel(Y)

            for i in range(x.size):
                s = data[:,y[i],x[i]]
                s = s / s.max()  # Normalize
                s = s - s.mean() # Remove mean to avoid triangular shape
                s = numpy.correlate(s, self.ref_s, mode='same')
                corr_cube[:,y[i],x[i]] = s

                temp = (((i + 1) * 100.00 / X.size))
                stdout.write('\r %2d%% ' % temp)
                stdout.write(self.loading[int(temp * 10 % 5)])
                stdout.flush()

            self.print(" Done.")
            corr_name = splitext(self.input_file)[0] + '--corrcube.fits'
            self.print(" Saving correlation cube to %s" % corr_name)

            corr_hdr = self.header.copy()
            corr_hdr.set('CORRFROM', self.input_file,'Cube used for corrcube.')
            corr_hdr.set('', '', before='CORRFROM')
            corr_hdr.set('', '--- Correlation cube ---', before='CORRFROM')

            writeto(corr_name, corr_cube, corr_hdr, clobber=True)
            del corr_hdr
            del corr_cube

            return corr_name

    def save(self):
        """
        Save files.
        """
        from astropy.io.fits import getdata, writeto
        from os.path import splitext

        v = self.verbose
        f = splitext(self.input_file)[0]
        h = self.header.copy()
        h.set('PHMREFX', self.ref_x)
        h.set('PHMREFY', self.ref_y)
        h.set('', '', before='PHMREFX')
        h.set('', '--- PHM Xtractor ---', before='PHMREFX')

        h.set('PHMREFF', self.input_file, 'Original file')
        h.set('PHMTYPE', 'observed')
        h.set('PHMUNIT', self.units)
        h.set('PHMSAMP', self.sample)

        filename = safe_save(f + "--obs_phmap.fits", overwrite=True, verbose=v)
        self.print(" Saving observed phase-map to file: %s" % filename)
        writeto(filename, self.phase_map, h, clobber=True)

        ## TODO Fix refspec file's header to keep calibration
        filename = safe_save(f + "--ref_spec.fits", overwrite=True, verbose=v)
        self.print(" Saving reference spectrum to file: %s" % filename)
        writeto(filename, self.ref_s, h, clobber=True)

        return

#==============================================================================
class PhaseMap_FP(PhaseMap):

    def __init__(self, filename, correlation=False, show=False, verbose=False):

        PhaseMap.__init__(self, filename, correlation=correlation,
                          show=show, verbose=verbose)

        # This is a Fabry-Perot data-cube. Let's make that clear to the user
        if self.verbose:
            print("\n Extracting phase-map from a Fabry-Perot data-cube.")

        # Measure the free-spectral-range
        self.free_spectral_range = self.get_free_spectral_range()

        # Getting reference spectrum
        self.ref_x, self.ref_y = self.find_reference_pixel()
        self.ref_s = self.get_reference_spectrum()

        # Calculate the FWHM
        self.fwhm = self.get_fwhm()

        # Calculate the finesse
        self.finesse = self.get_finesse()

        if self.verbose:
            print(" Ideal number of channels: %.1f channels"
                  % round(2 * self.finesse))
            print(" Ideal sampling: %.1f %s / channel"
                  % (self.free_spectral_range / round(2 * self.finesse),
                     self.units))

        if correlation:
            self.extract_from = self.use_correlation()
        else:
            self.extract_from = self.input_file

        self.phase_map = self.extract_phase_map()
        self.save()

        return

    def extract_phase_map(self):
        """
        Extract the phase-map.
        """

        now = time.time()

        try:
            sampling = self.header['C3_3']
        except KeyError:
            sampling = 1

        fsr = round(self.free_spectral_range / sampling)

        # Reading data
        if self.verbose:
            print("\n Starting phase-map extraction.")
            print(" Reading data from %s file" % self.extract_from)
        data = pyfits.getdata(self.extract_from)
        # data = data[0:fsr]

        # Extracting phase-map
        if self.verbose:
            print(" Extracting phase-map...")
        data = numpy.where(data > data.mean() + data.std(), data, -numpy.inf)
        phase_map = numpy.argmax(data, axis=0) * sampling

        if self.verbose:
            print(" Done in %.2f seconds" % (time.time() - now))
        return phase_map

    def find_reference_pixel(self):
        """
        Read the reference pixel from header or find it.
        """
        if self.verbose:
            print("\n Finding reference pixel.")

        if ('PHMREFX' in self.header) and ('PHMREFY' in self.header):

            if self.verbose:
                print(" Found reference pixel found in header.")

            ref_x = self.header['PHMREFX']
            ref_y = self.header['PHMREFY']

            if self.verbose:
                print(" Using [%d, %d]" % (ref_x, ref_y))

        else:
            if self.verbose:
                print(" Reference pixel NOT found in header.")
                print(" Trying to find the center of the rings.")
            ref_x, ref_y = self.find_rings_center()

        return ref_x, ref_y

    def find_rings_center(self):
        """
        Method used to find the center of the rings inside a FP data-cube.
        """

        now = time.time()

        # Renaming some variables
        width = self.width
        height = self.height
        fsr = round(self.free_spectral_range / self.header['C3_3'])

        # Choosing the points
        x = (numpy.linspace(0.2, 0.8, 500) * width).astype(int)
        y = (numpy.linspace(0.2, 0.8, 500) * height).astype(int)

        ref_x = self.header['NAXIS1'] // 2
        ref_y = self.header['NAXIS2'] // 2

        self.print(" Start center finding.")
        old_ref_x = ref_x
        old_ref_y = ref_y

        if self.show:
            pyplot.figure()

        for i in range(6):

            ref_y = max(ref_y, 0)
            ref_y = min(ref_y, self.header['NAXIS2'])

            ref_x = max(ref_x, 0)
            ref_x = min(ref_x, self.header['NAXIS2'])

            temp_x = self.data[:fsr, ref_y, x]
            temp_y = self.data[:fsr, y, ref_x]

            temp_x = numpy.argmax(temp_x, axis=0)
            temp_y = numpy.argmax(temp_y, axis=0)

            px = scipy.polyfit(x, temp_x, 2)
            py = scipy.polyfit(y, temp_y, 2)

            ref_x = round(- px[1] / (2.0 * px[0]))
            ref_y = round(- py[1] / (2.0 * py[0]))

            if self.show:
                pyplot.title("Finding center of the rings")
                pyplot.cla()
                pyplot.plot(x, temp_x, 'b.', alpha=0.25)
                pyplot.plot(x, scipy.polyval(px, x), 'b-', lw=2)
                pyplot.plot(y, temp_y, 'r.', alpha=0.25)
                pyplot.plot(y, scipy.polyval(py, y), 'r-', lw=2)
                pyplot.gca().yaxis.set_ticklabels([])
                pyplot.axvline(ref_x, ls='--', c='blue', label='x')
                pyplot.axvline(ref_y, ls='--', c='red', label='y')
                pyplot.legend(loc='best')
                pyplot.grid()
                pyplot.ylabel("Iteration number %d" %(i+1))

            # Selecting valid data
            error_x = numpy.abs(temp_x - scipy.polyval(px, x))
            error_y = numpy.abs(temp_y - scipy.polyval(py, y))

            cond_x = numpy.where(error_x <= 3 * error_x.std(), True, False)
            cond_y = numpy.where(error_y <= 3 * error_y.std(), True, False)

            x = x[cond_x]
            y = y[cond_y]

            # Choosing when to stop
            if (abs(old_ref_x - ref_x) <= 2) and (abs(old_ref_y - ref_y) <= 2):

                try:
                    # If the cube was binned this will be useful
                    ref_x = (ref_x - self.header['CRPIX1'] + 1) \
                            * self.header['CDELT1'] + self.header['CRVAL1']

                    # If the cube was binned this will be useful
                    ref_y = (ref_y - self.header['CRPIX2']) \
                            * self.header['CDELT2'] + self.header['CRVAL2']
                except KeyError:
                    pass

                if self.verbose:
                    print(" Rings center found at: [%d, %d]" % (ref_x, ref_y))
                    print(" Done in %.2f s" % (time.time() - now))

                if self.show:
                    pyplot.tight_layout()
                    pyplot.show()

                return ref_x, ref_y

            else:
                old_ref_x = ref_x
                old_ref_y = ref_y

        if self.show:
            pyplot.tight_layout()
            pyplot.show()

        if self.verbose:
            print(" Rings center NOT found.")

        ref_x = self.header['NAXIS1'] // 2
        ref_y = self.header['NAXIS2'] // 2

        # If the cube was binned this will be useful
        try:
            ref_x = (ref_x - self.header['CRPIX1']) \
                    * self.header['CDELT1'] + self.header['CRVAL1']
            ref_y = (ref_y - self.header['CRPIX2']) \
                    * self.header['CDELT2'] + self.header['CRVAL2']
        except:
            pass

        if self.verbose:
            print(" Done in %.2f s" % (time.time() - now))
            print(" Using [%d, %d]." % (ref_x, ref_y))
        return ref_x, ref_y

    def get_finesse(self):
        """
        Assuming you have the Free-Spectral-Range in Z unit and that
        you have the FWHM in Z units as well, calculate the finesse.
        """

        finesse = self.free_spectral_range / self.fwhm

        if self.verbose:
            print(" Finesse = %.1f" % finesse)

        return finesse


    def get_free_spectral_range(self):
        """
        A quick-and-dirty way to measure the free range in FP units.
        The method subtracts each frame of the data-cube from the
        first one. Then, it calculates the absolute value and collapse
        in X and Y. The FSR is where the resulting spectrum is minimum,
        excluding (of course), the first one.
        """

        if self.verbose:
            print(" Finding the free-spectral-range.")

        now = time.time()

        # First frame is the reference frame
        ref_frame = self.data[0,:,:]

        # Subtract all frames from the first frame
        data = self.data - ref_frame

        # Get the absolute value
        data = numpy.abs(data)

        # Sum over the spatial directions
        data = data.sum(axis=2)
        data = data.sum(axis=1)

        # Interpolate data
        s = interpolate.UnivariateSpline(self.z, data, k=3)
        z = numpy.linspace(self.z[5:].min(), self.z.max(), 1000)

        # Find the free-spectral-range in z units
        fsr = z[numpy.argmin(s(z))] - self.z[0]

        # Find the free-spectral-range in number of channels
        fsr_channel = numpy.argmin(numpy.abs(self.z - z[numpy.argmin(s(z))]))

        # Calculate the sampling
        sampling = fsr / fsr_channel

        if self.verbose:
            print(" FSR = %.1f %s" % (fsr, self.units))
            print("     = %d channels" % fsr_channel)
            print(" Sampling = %.1f %s / channel" % (sampling, self.units))
            print(" Done in %.2f s" % (time.time() - now))

        # Plot to see how it goes
        if self.show:
            pyplot.title("Finding the Free-Spectral-Range")
            pyplot.plot(self.z, data, 'bo', label='Measured data')
            pyplot.plot(z, s(z), 'r-', lw=2, label='3rd deg spline fitting')
            pyplot.xlabel("z [%s]" % self.units)
            pyplot.axvline(x=(fsr + self.z.min()), ls='--', c='red',
                           label='Free-Spectral-Range \nat z = %.1f' % fsr)
            pyplot.legend(loc='best')
            pyplot.gca().yaxis.set_ticklabels([])
            pyplot.grid()
            pyplot.tight_layout()
            pyplot.show()

        return fsr

    def save(self):
        """
        Save files.
        """
        from os.path import splitext

        v = self.verbose
        f = splitext(self.input_file)[0]
        h = self.header.copy()
        h['PHMREFX'] = self.ref_x
        h['PHMREFY'] = self.ref_y
        h.set('', '', before='PHMREFX')
        h.set('', '--- PHM Xtractor ---', before='PHMREFX')

        # if 'PHMREFX' not in self.header:
        #     update = '.'
        #     while update.upper() not in 'YESNO':
        #         update = raw_input(" Update input file? [Y]/n \n ")
        #         if update.upper() in 'YES':
        #             self.print(" Updating input file %s" % self.input_file)
        #             data = getdata(self.input_file)
        #             writeto(self.input_file, data, h, clobber=True)

        fsr = self.free_spectral_range

        h['PHMREFF'] = (self.input_file, 'Original file')
        h['PHMTYPE'] = 'observed'
        h['PHMUNIT'] = self.header['CUNIT3']
        h['PHMFSR'] = (round(fsr, 2),
                       'Free-spectral-range in %s units' % self.units)
        h['PHMSAMP'] = (self.header['C3_3'], 'Used sample [%s / channel].'
                                             % self.units)

        # TODO Remove 3rd axis calibration residuals
        filename = safe_save(f + "--obs_phmap.fits", overwrite=True, verbose=v)
        self.print(" Saving observed phase-map to file: %s" % filename)
        pyfits.writeto(filename, self.phase_map, h, clobber=True)

        filename = safe_save(f + "--ref_spec.fits", overwrite=True, verbose=v)
        self.print(" Saving reference spectrum to file: %s" % filename)
        pyfits.writeto(filename, self.ref_s, h, clobber=True)

        return

#==============================================================================
class PhaseMap_iBTF(PhaseMap):

    def __init__(self, filename, correlation=False, show=False, verbose=False):

        PhaseMap.__init__(self, filename, correlation=correlation,
                          show=show, verbose=verbose)

        # This is an iBTF data-cube. Let's make that clear to the user
        if self.verbose:
            print("\n  Extracting phase-map from a iBTF data-cube.")

        # Getting reference spectrum
        self.ref_x, self.ref_y = self.find_reference_pixel()
        self.ref_s = self.get_reference_spectrum()

        # Use correlation?
        if correlation:
            self.extract_from = self.use_correlation()
        else:
            self.extract_from = self.input_file

        # Extract phase-map
        self.phase_map = self.extract_phase_map()
        self.save()

#==============================================================================
if __name__ == '__main__':
    main()