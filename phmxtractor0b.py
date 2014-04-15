#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    v01a - Original File.
    v01b - Added '-s'/'--show' option.
         - "--correlation" option default set to FALSE.
         - Changed all the astropy.io.fits imports to a single import at the
           beginning of the file.
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
    if v:
        start = time.time()
        print("\n Phase-Map Extractor")
        print(" by Bruno Quint & Fabricio Ferrari")
        print(" version 0.1b - Apr 2014")
        print("\n Extracting phase-map from file: %s" % args.filename)

    # Checking input data -----------------------------------------------------
    # TODO Add a manual option for the case where the instrument was not recognized.
    if v: 
        print(" Checking data-cube for phase-correction.")
        
    if not is_btfi_data(args.filename):
        print(" Instrument type not recognized.")
        print(" Leaving now.")
        sys.exit()
        
    # Selecting BTFI mode and extracting phase-map -----------------------------
    mode = get_btfi_mode(args.filename)
    if mode == 'ibtf':
        PhaseMap_iBTF(args.filename, correlation=args.correlation, verbose=v)
    elif mode == 'fp':
        PhaseMap_FP(args.filename, correlation=args.correlation,
                    show=args.show, verbose=v)
    
    # All done! ---------------------------------------------------------------
    if v:
        print(" Done.")
        end = time.time() - start
        print("\n Total time ellapsed: %02d:%02d:%02d" % 
              (end // 3600, end % 3600 // 60, end % 60))     
        print(" All done!\n")

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

def get_btfi_mode(filename):
    """
    Return if BTFI was obtained with a Fabry-Perot or with the iBTF.
    """
    from astropy.io.fits import getheader
    
    header = getheader(filename)
    
    if header['INSTRMOD'].upper() in ['IBTF']:
        return 'ibtf'
    
    if header['INSTRMOD'].upper() in ['FP', 'FABRY-PEROT']:
        return 'fp'
    
def is_btfi_data(filename):
    """
    Check if input file was obtained with BTFI instrument.
    """
    from astropy.io.fits import getheader
    
    header = getheader(filename)
    btfi_data = ('INSTRUME' in header)
    btfi_data = btfi_data and (header['INSTRUME'].upper() in ['BTFI'])
    return btfi_data  


#==============================================================================
class PhaseMap:

    def __init__(self, filename, correlation=False, show=False, verbose=False):

        # Setting main configuration ------------------------------------------
        self.input_file = filename
        self.verbose = verbose
        self.show = show

        # Reading raw data ----------------------------------------------------
        self.print(" Loading data.")
        self.data = pyfits.getdata(filename)
        self.header = pyfits.getheader(filename)
        self.print(" Done.")

        # Reading data-cube configuration -------------------------------------
        self.width = self.header['NAXIS1']
        self.height = self.header['NAXIS2']
        self.depth = self.header['NAXIS3']

        # Reading Z calibration for plotting ----------------------------------
        self.z = self.get_calibration()
        self.units = self.header['CUNIT3']

        return
    
    def extract_phase_map(self):
        """
        Extract the phase-map.
        """
        from astropy.io.fits import getdata
        from numpy import argmax, inf, where
        
        self.print("\n Starting phase-map extraction.")
        self.print(" Reading data from %s file" % self.extract_from)
        data = getdata(self.extract_from)
        data = where(data > data.mean() + data.std(), data, -inf)
        phase_map = argmax(data, axis=0) * self.header['CDELT3']
        
        return phase_map
    
    def find_reference_pixel(self):
        """Read the reference pixel from header or find it."""
        if ('PHMREFX' in self.header) and ('PHMREFY' in self.header):
            self.print(" \n Found reference pixel in header.")
            ref_x = self.header['PHMREFX']
            ref_y = self.header['PHMREFY']
            self.print(" Using [%d, %d]" % (self.ref_x, self.ref_y))
        else:
            self.print(" \n Reference pixel NOT in header.")
            ref_x = self.width // 2
            ref_y = self.height // 2
            self.print(" Using [%d, %d]" % (self.ref_x, self.ref_y))
            
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
            print("[!] Calibration in third axis not found.")
            print("[!] I will ignore this step.")

        return z


    def get_reference_spectrum(self):
        """
        Get the reference spectrum.
        """ 
        ref_s = pyfits.getdata(self.input_file)[:,self.ref_y, self.ref_x]
        ref_s = ref_s / ref_s.max()  # Normalize
        ref_s = ref_s - ref_s.mean() # Remove mean to avoid triangular shape

        if self.show:
            pyplot.figure()
            pyplot.title("Reference Spectrum")
            pyplot.plot(self.z, ref_s, 'ko-', label="Reference spectrum")
            pyplot.grid()
            pyplot.xlabel("z [%s]" % self.units)
            pyplot.gca().set_yticklabels([])
            pyplot.show()
        
        return ref_s
            
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
        
        if 'PHMREFX' not in self.header:
            update = '.'
            while update.upper() not in 'YESNO':
                update = raw_input(" Update input file? [Y]/n \n ")
                if update.upper() in 'YES':
                    self.print(" Updating input file %s" % self.input_file)
                    data = getdata(self.input_file)
                    writeto(self.input_file, data, h, clobber=True)        
                    
        h.set('PHMREFF', self.input_file, 'Original file')
        h.set('PHMTYPE', 'observed')
        h.set('PHMUNIT', self.header['CUNIT3'])
        
        filename = safe_save(f + "--obs_phmap.fits", overwrite=True, verbose=v)
        self.print(" Saving observed phase-map to file: %s" % filename)
        writeto(filename, self.phase_map, h, clobber=True)
        
        filename = safe_save(f + "--ref_spec.fits", overwrite=True, verbose=v)
        self.print(" Saving reference spectrum to file: %s" % filename)
        writeto(filename, self.ref_s, h, clobber=True)
        
        return
                
#==============================================================================
class PhaseMap_FP(PhaseMap):

    def __init__(self, filename, correlation=False, show=False, verbose=False):

        PhaseMap.__init__(self, filename, correlation=correlation,
                          show=show, verbose=verbose)

        # This is a Fabry-Perot data-cube. Let's make that clear to the user --
        self.print("\n Fabry-Perot data ---------------")

        # Measure the free-spectral-range -------------------------------------
        self.free_spectral_range = self.get_free_spectral_range()

        # Getting reference spectrum ------------------------------------------
        self.ref_x, self.ref_y = self.find_reference_pixel()
        self.ref_s = self.get_reference_spectrum()

        # Calculate the FWHM --------------------------------------------------
        self.fwhm = self.get_fwhm()

        # if correlation:
        #     self.extract_from = self.use_correlation()
        # else:
        #     self.extract_from = self.input_file
        #
        # self.phase_map = self.extract_phase_map()
        # self.save()


    def find_reference_pixel(self):
        """
        Read the reference pixel from header or find it.
        """

        if ('PHMREFX' in self.header) and ('PHMREFY' in self.header):
            self.print(" \n Found reference pixel found in header.")
            ref_x = self.header['PHMREFX']
            ref_y = self.header['PHMREFY']
            self.print(" Using [%d, %d]" % (ref_x, ref_y))

        else:
            self.print(" \n Reference pixel NOT found in header.")
            self.print(" Trying to find the center of the rings.")
            ref_x, ref_y = self.find_rings_center()

        return ref_x, ref_y


    def find_rings_center(self):
        """
        Method used to find the center of the rings inside a FP data-cube.
        """

        # Renaming some variables ---------------------------------------------
        width = self.width
        height = self.height
        fsr = round(self.free_spectral_range / self.header['C3_3'])
        
        # Choosing the points -------------------------------------------------
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
                pyplot.subplot(3,1,i+1)
                pyplot.plot(x, temp_x, 'b.', alpha=0.25)
                pyplot.plot(x, scipy.polyval(px, x), 'b-', lw=2)
                pyplot.plot(y, temp_y, 'r.', alpha=0.25)
                pyplot.plot(y, scipy.polyval(py, y), 'r-', lw=2)
                pyplot.gca().yaxis.set_ticklabels([])
                pyplot.axvline(ref_x, ls='--', c='blue', label='x')
                pyplot.axvline(ref_y, ls='--', c='red', label='y')
                pyplot.legend(loc='best')
                pyplot.grid()
                pyplot.xlabel("Iteration number %d" %(i+1))
            
            # Selecting valid data
            error_x = numpy.abs(temp_x - scipy.polyval(px, x))
            error_y = numpy.abs(temp_y - scipy.polyval(py, y))
                        
            cond_x = numpy.where(error_x <= 3 * error_x.std(), True, False)
            cond_y = numpy.where(error_y <= 3 * error_y.std(), True, False)
            
            x = x[cond_x]
            y = y[cond_y]
            
            # Choosing when to stop
            if (abs(old_ref_x - ref_x) <= 2) and (abs(old_ref_y - ref_y) <= 2):
                
                ref_x = (ref_x - self.header['CRPIX1'] + 1) \
                        * self.header['CDELT1'] + self.header['CRVAL1']

                ref_y = (ref_y - self.header['CRPIX2']) \
                        * self.header['CDELT2'] + self.header['CRVAL2']

                if self.verbose:
                    print(" Rings center found at: [%d, %d]" % (ref_x, ref_y))

                if self.show:
                    pyplot.show()

                return ref_x, ref_y
            
            else:
                old_ref_x = ref_x
                old_ref_y = ref_y

        if self.show:
            pyplot.show()
        print(" Rings center NOT found.")
        
        ref_x = self.header['NAXIS1'] // 2
        ref_y = self.header['NAXIS2'] // 2
        
        ref_x = (ref_x - self.header['CRPIX1']) * self.header['CDELT1'] + self.header['CRVAL1']
        ref_y = (ref_y - self.header['CRPIX2']) * self.header['CDELT2'] + self.header['CRVAL2']
        
        print(" Using [%d, %d]." % (ref_x, ref_y))        
        return ref_x, ref_y


    def get_free_spectral_range(self):
        """
        Method created to find the Free-Spectral-Range of a calibration
        data-cube. Be sure that you scanned more than a single FRS.
        """
        now = time.time()
        self.print(" Finding free-spectral-range...")

        if self.show:
            pyplot.figure()

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

        # Find the free-spectral-range
        fsr = z[numpy.argmin(s(z))] - self.z.min()

        if self.verbose:
            print("  FSR = %.1f %s" % (fsr, self.units))
            print("  Done in %.2f s" % (time.time() - now))

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
            pyplot.show()

        return fsr

    def get_fwhm(self):
        """
        Returns the full-width-at-half-maximum.
        """
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
        fwhm1 = 2.35482 * p[2]

        zzz = zz[sss > 0]
        fwhm2 = zzz.ptp()

        if self.show:
            pyplot.figure()
            pyplot.title("Measure the FWHM")
            pyplot.plot(z, s, 'bo')
            pyplot.plot(zz, ss(zz), 'b-', lw=2)
            pyplot.plot(zz, sss, 'r-', lw=2, alpha=0.3)
            pyplot.plot(zz, fit_func(p, zz), 'g-', lw=2, alpha=0.3)
            pyplot.axvline(p[1] - fwhm1 / 2, ls='--', c='green', lw=2)
            pyplot.axvline(p[1] + fwhm1 / 2, ls='--', c='green', lw=2,
                label='Gauss Fit = %.1f %s' % (fwhm1, self.units))
            pyplot.axvline(p[1] + fwhm2 / 2, ls='--', c='red', lw=2)
            pyplot.axvline(p[1] - fwhm2 / 2, ls='--', c='red', lw=2,
                label='Definition = %.1f %s' % (fwhm2, self.units))
            pyplot.legend(loc='best')
            pyplot.grid()
            pyplot.show()

        return fwhm2


#==============================================================================
class PhaseMap_iBTF(PhaseMap):
    pass
        
#==============================================================================
if __name__ == '__main__':
    main()