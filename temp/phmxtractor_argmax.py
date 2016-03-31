#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import division, print_function

import argparse
import matplotlib.pyplot as pyplot
import multiprocessing
import numpy
import pyfits
import time
import scipy.interpolate
import scipy.optimize
import scipy.ndimage
import sys

def main():
    
    # Parse arguments ---------------------------------------------------------
    parser = argparse.ArgumentParser(description="Extracts the phase-map" + 
                                     "from a fits file containing a data-cube.")
    
    parser.add_argument('-c', '--correlation', action='store_true',
                        help="Use correlation cube? true/[false]")
    
    parser.add_argument('filename', type=str, help="Input data-cube name.")
    
    parser.add_argument('-o', '--output', type=str, default=None,
                        help="Name of the output phase-map file.")
    
    parser.add_argument('-q', '--quiet', action='store_true', 
                        help="Run program quietly.")

    parser.add_argument('-s', '--show_plots', action='store_true',
                        help="Show plots involved in the process.")
    
    args = parser.parse_args()
    
    # Starting program --------------------------------------------------------
    v = not args.quiet
    if v:
        start = time.time()
        print("\n Phase-Map Extractor")
        print(" by Bruno Quint & Fabricio Ferrari")
        print(" version 0.1a - Jan 2014")
        print("\n Extracting phase-map from file: %s" % args.filename)

    # Checking input data -----------------------------------------------------
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
                    show=args.show_plots, verbose=v)
    
    # All done! ---------------------------------------------------------------
    if v:
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
    from pyfits import getheader
    
    header = getheader(filename)
    
    if header['INSTRMOD'].upper() in ['IBTF']:
        return 'ibtf'
    
    if header['INSTRMOD'].upper() in ['FP', 'FABRY-PEROT']:
        return 'fp'
    
def is_btfi_data(filename):
    """
    Check if input file was obtained with BTFI instrument.
    """
    from pyfits import getheader
    
    header = getheader(filename)
    btfi_data = ('INSTRUME' in header)
    btfi_data = btfi_data and (header['INSTRUME'].upper() in ['BTFI'])
    return btfi_data  

class PhaseMap:
    def __init__(self, filename, correlation=False, show=False, verbose=False):
        """
        Creating phase-map.
        """
        self.input_file = filename
        self.data = pyfits.getdata(filename)
        self.header = pyfits.getheader(filename)
        self.show = show
        self.verbose = verbose
        self.loading = [' ','-','\\','|','/']
        
        self.width = self.header['NAXIS1']
        self.height = self.header['NAXIS2']
        self.depth = self.header['NAXIS3']
        self.units = self.header['CUNIT3']

        self.zref = self.header['CRPIX3'] - 1
        self.zval = self.header['CRVAL3']
        self.zdelta = self.header['CDELT3']
        
        self.print(" File obtained through an %s scan." \
                   % self.header['INSTRMOD'])

        return
    
    def extract_phase_map(self):
        """
        Extract the phase-map.
        """
        from functools import partial

        now = time.time()
        self.print("\n Starting phase-map extraction.")
        self.print(" Reading data from %s file" % self.input_file)

        z = numpy.argmax(self.data, axis=0)
        z = z - self.zref
        z = z * self.zdelta
        z = z + self.zval
        phase_map = z

        now = time.time() - now
        self.print(" Time spent: %02dh%02dm%02.2fs" %
                      (now // 3600, now % 3600 // 60, now % 60))
        print(" Done.")
        return phase_map

    def find_reference_pixel(self):
        """
        Read the reference pixel from header or find it.
        """
        if ('PHMREFX' in self.header) and ('PHMREFY' in self.header):
            self.print(" \n Found reference pixel in header.")
            ref_x = self.header['PHMREFX']
            ref_y = self.header['PHMREFY']
            self.print(" Using [%d, %d]" % (ref_x, ref_y))
        else:
            self.print(" \n Reference pixel NOT in header.")
            ref_x = self.width // 2
            ref_y = self.height // 2
            self.print(" Using [%d, %d]" % (ref_x, ref_y))
            
        return ref_x, ref_y
    
    def get_reference_spectrum(self):
        """
        Get the reference spectrum.
        """
        ref_s = self.data[:,self.ref_y, self.ref_x]
        ref_s = ref_s / ref_s.max()  # Normalize
        ref_s = ref_s - ref_s.mean() # Remove mean to avoid triangular shape

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
        
        from pyfits import getdata, getheader, writeto
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
        from pyfits import getdata, writeto
        from os.path import splitext
        
        v = self.verbose
        f = splitext(self.input_file)[0]
        h = self.header.copy()
        h.set('PHMREFX', self.ref_x)
        h.set('PHMREFY', self.ref_y)
        h.set('', '', before='PHMREFX')
        h.set('', '--- PHM Xtractor ---', before='PHMREFX')
        
        #if 'PHMREFX' not in self.header:
        #    update = '.'
        #    while update.upper() not in 'YESNO':
        #        update = raw_input(" Update input file? [Y]/n \n ")
        #        if update.upper() in 'YES':
        #            self.print(" Updating input file %s" % self.input_file)
        #            data = getdata(self.input_file)
        #            writeto(self.input_file, data, h, clobber=True)        
                    
        h.set('PHMREFF', self.input_file, 'Original file')
        h.set('PHMTYPE', 'observed')
        h.set('PHMUNIT', self.header['CUNIT3'])
        h.set('PHMSAMP', self.header['CDELT3'])
        
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

        self.free_range_in_z = self.find_free_range()
        self.free_range_in_channels = self.free_range_in_z / self.zdelta
        self.data = self.data[0:self.free_range_in_channels+1,:,:]
        self.header["PHMFSR"] = (self.free_range_in_z, "Free range in %s" % self.units)

        self.ref_x, self.ref_y = self.find_reference_pixel()
        self.ref_s = self.get_reference_spectrum()
        self.z_fwhm = self.get_fwhm()
        # self.data = self.data[0:int(self.z_fwhm / self.zdelta),:,:]

        # if correlation:
        #     self.extract_from = self.use_correlation()
        # else:
        #     self.extract_from = self.input_file

        self.phase_map = self.extract_phase_map()
        self.save()

    def find_free_range(self):
        """
        A quick-and-dirty way to measure the free range in FP units.
        The method subtracts each frame of the data-cube from the
        first one. Then, it calculates the absolute value and colapse
        in X and Y. The FSR is where the resulting spectrum is minimum,
        excluding (of course), the first one.
        """

        if self.verbose:
            print(" Finding the free-range.")

        now = time.time()
        ref_frame = self.data[0].reshape((1, self.height, self.width))
        s = self.data - numpy.repeat(ref_frame, self.depth, axis=0)
        s = numpy.abs(s)
        s = numpy.average(s, axis=2)
        s = numpy.average(s, axis=1)

        z = (numpy.arange(self.depth) - self.zref) * self.zdelta + self.zval
        f = scipy.interpolate.interp1d(z, s, kind='cubic')

        zz = numpy.linspace(z[5], z[-1], 1000)
        ss = f(zz)

        min_z = zz[numpy.argmin(ss)]
        min_channel = numpy.argmin(numpy.abs(z - min_z))
        free_range = min_z - z[0] + self.zdelta
        used_sample = free_range / min_channel

        now = time.time() - now
        if self.verbose:
            print(" FSR = %d %s / %d channels" % (free_range, self.units, min_z+1))
            print(" Sampling used = %.1f %s / channel" % (used_sample, self.units))
            print(" Time spent: %02dh%02dm%02.2fs" %
                  (now // 3600, now % 3600 // 60, now % 60))
            print(" Done.")

        if self.show:
            pyplot.figure()
            pyplot.plot(z, s, 'bo')
            pyplot.plot(zz, ss, 'r--')
            pyplot.axvline(zz[numpy.argmin(ss)], color='red')
            pyplot.axvline(z[0], color='red')
            pyplot.gca().yaxis.set_ticklabels([])
            pyplot.xlabel("z [%s]" % self.units)
            pyplot.grid()
            pyplot.title("Find the Free-Range")
            pyplot.show()

        self.header['PHM_FSR'] = (free_range,
            'Free-spectral-range in %s units' % self.units)
        self.header['PHM_SMP'] = (used_sample,
                                  'Used sample [%s / channel].' % self.units)

        return free_range

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
        # Renaming some variables
        width = self.width
        height = self.height

        # Choosing the points
        x = (numpy.linspace( 0.2, 0.8, 500) * width).astype(int)
        y = (numpy.linspace( 0.2, 0.8, 500) * height).astype(int)

        ref_x = self.header['NAXIS1'] // 2
        ref_y = self.header['NAXIS2'] // 2
        data = self.data

        self.print(" Start center finding.")
        now = time.time()
        old_ref_x = ref_x
        old_ref_y = ref_y

        if self.show:
            pyplot.figure(figsize=(12,4))

        for i in range(10):

            temp_x = data[:, ref_y, x]
            temp_y = data[:, y, ref_x]

            temp_x = numpy.argmax(temp_x, axis=0)
            temp_y = numpy.argmax(temp_y, axis=0)

            # Fitting parabola
            px = scipy.polyfit(x, temp_x, 2)
            py = scipy.polyfit(y, temp_y, 2)

            # Peak finding using parabola equation
            ref_x = round(- px[1] / (2.0 * px[0]))
            ref_y = round(- py[1] / (2.0 * py[0]))

            if self.show:
                pyplot.subplot(1,3,(i+1) % 3)
                if i == 1:
                    pyplot.title("Finding rings center")
                pyplot.plot(x, temp_x, 'b.', label='x')
                pyplot.plot(x, scipy.polyval(px, x), 'b-')
                pyplot.plot(y, temp_y, 'r.', label='y')
                pyplot.plot(y, scipy.polyval(py, y), 'r-')
                pyplot.legend(loc='best')
                pyplot.gca().xaxis.set_ticklabels([])
                pyplot.gca().yaxis.set_ticklabels([])
                pyplot.grid()
                pyplot.xlabel("Interaction %d" % (i + 1))

            # Selecting valid data
            error_x = numpy.abs(temp_x - scipy.polyval(px, x))
            error_y = numpy.abs(temp_y - scipy.polyval(py, y))

            cond_x = numpy.where(error_x <= error_x.std(), True, False)
            cond_y = numpy.where(error_y <= error_y.std(), True, False)

            x = x[cond_x]
            y = y[cond_y]

            # Choosing when to stop -------------------------------------------
            if (abs(old_ref_x - ref_x) <= 2) and (abs(old_ref_y - ref_y) <= 2):

                try:
                    ref_x = (ref_x - self.header['CRPIX1']) * self.header['CDELT1'] + self.header['CRVAL1']
                    ref_y = (ref_y - self.header['CRPIX2']) * self.header['CDELT2'] + self.header['CRVAL2']
                except KeyError:
                    self.print(" [!] No spatial calibration found. Using pixel values.")

                now = time.time() - now
                self.print(" Rings center found at: [%d, %d]" % (ref_x, ref_y))
                self.print(" Time spent: %02dh%02dm%02.2fs" %
                      (now // 3600, now % 3600 // 60, now % 60))
                print(" Done.")

                if self.show:
                    pyplot.show()

                return ref_x, ref_y

            else:
                old_ref_x = ref_x
                old_ref_y = ref_y

        self.print(" Rings center NOT found.")

        ref_x = self.header['NAXIS1'] // 2
        ref_y = self.header['NAXIS2'] // 2

        ref_x = (ref_x - self.header['CRPIX1']) * self.header['CDELT1'] + self.header['CRVAL1']
        ref_y = (ref_y - self.header['CRPIX2']) * self.header['CDELT2'] + self.header['CRVAL2']

        now = time.time() - now
        self.print(" Time used for that %02dh%02dm%02.2fs" %
                   (now // 3600, now % 3600 // 60, now % 60))
        self.print(" Using [%d, %d]." % (ref_x, ref_y))
        return ref_x, ref_y

    def get_fwhm(self):

        now = time.time()
        self.print("\n Measuring the FWHM in z in %s units" % self.units)
        z = numpy.arange(self.data.shape[0]) - self.zref
        z = z * self.zdelta
        z = z + self.zval
        s = self.ref_s
        #s = scipy.ndimage.median_filter(s, 3)

        avg = s.mean()
        sig = s.std()
        cond = s.max() / 2.0

        mean = (z * s)[s > cond].sum() / s[s > cond].sum()
        fwhm = numpy.sqrt((s * (z - mean) ** 2)[s > cond].sum()
                          / s[s > cond].sum())

        fit_func = lambda p, x: p[0] * numpy.exp(-(p[1] - x) ** 2 / (2 * p[2] ** 2))
        err_func = lambda p, x, y: y - fit_func(p, x)
        p, _ = scipy.optimize.leastsq(err_func, [s.max(), mean, fwhm], (z, s))
        zz = numpy.linspace(z[0], z[-1], 1000)

        z_central = p[1]
        fwhm = p[2] * (2. * numpy.sqrt(2. * numpy.log(2)))
        finesse = (self.free_range_in_z / fwhm)
        n_channels = 2 * finesse
        steps = self.free_range_in_z / n_channels
        now = time.time() - now

        self.print(" The fwhm is: %.1f %s" % (fwhm, self.units))
        self.print(" The finesse is: %.1f" % finesse)
        self.print(" The minimum number of frames in a " +
                   "free-spectral-range is: %d" % round(n_channels))
        self.print(" The scanning step dZ shall be smaller than: %.1f %s" %
                   (steps, self.units))
        self.print(" Time spend: %02dh%02dm%02.2fs" % (
            now // 3600, now % 3600 // 60, now % 60 ))
        self.print(" Done.")

        if self.show:
            pyplot.figure()
            pyplot.plot(z, s, 'bo-')
            pyplot.plot(zz, fit_func(p, zz), 'r-')
            pyplot.grid()
            pyplot.xlabel("z [%s]" % self.units)
            pyplot.gca().yaxis.set_ticklabels([])
            pyplot.axvline(z_central, color='k', linestyle='-')
            pyplot.axvline(z_central - fwhm // 2, color='k', linestyle='--')
            pyplot.axvline(z_central + fwhm // 2, color='k', linestyle='--')
            pyplot.show()

        self.header['PHMFWHM'] = (fwhm, 'FWHM in %s units' % self.units)
        self.header['PHM_FIN'] = (finesse, 'Calculated finesse')

        return

                
#==============================================================================
class PhaseMap_iBTF(PhaseMap):
    pass
        
#==============================================================================
if __name__ == '__main__':
    main()
