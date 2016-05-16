#!/usr/bin/env python2
# -*- coding: utf8 -*-
"""
    SAMI XJoin

    This script simply joins the four existing extensions inside a FITS file
    created during observations with SAMI (SAM Imager). During the process,
    it also fits a 2nd degree polynomium to the OVERSCAN region that is
    subtracted from the corresponding image.

    The user also may want to add flags in order to process the images
    according to the following options (in order):

    - BIAS subtraction;
    - DARK subtraction;
    - Remove hot pixels and cosmic rays;
    - Remove overglow using a long exposure DARK image;
    - Divide by the FLAT;
    - Divide by the exposure time;

    The documentation for each process is shown in the corresponding function.

    Bruno Quint (bquint at ctio.noao.edu)
    May 2016

    Thanks to Andrei Tokovinin and Claudia M. de Oliveira for the ideas that
    were implemented here.
"""

from __future__ import division, print_function

import astropy.io.fits as pyfits
import argparse
import logging as log
import numpy as np
import os

from scipy import ndimage
from scipy import signal

try:
    xrange
except NameError:
    # noinspection PyShadowingBuiltins
    xrange = range

# Piece of code from cosmics.py
# We define the laplacian kernel to be used
laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])

# Other kernels :
growkernel = np.ones((3, 3))

# dilation structure for some morphological operations
dilstruct = np.ones((5, 5))
dilstruct[0, 0] = 0
dilstruct[0, 4] = 0
dilstruct[4, 0] = 0
dilstruct[4, 4] = 0


# noinspection PyPep8Naming
class SAMI_XJoin:
    def __init__(self, list_of_files, bias_file=None, clean=False,
                 cosmic_rays=False, dark_file=None, debug=False,
                 flat_file=None, glow_file=None, time=False, verbose=False):

        self.set_verbose(verbose)
        self.set_debug(debug)
        self.main(list_of_files, bias_file=bias_file, clean=clean,
                  cosmic_rays=cosmic_rays, dark_file=dark_file,
                  flat_file=flat_file, glow_file=glow_file, time=time)

        return

    @staticmethod
    def bias_subtraction(data, header, prefix, bias_file):

        if bias_file is not None:
            bias = pyfits.getdata(bias_file)
            data -= bias
            header['BIASFILE'] = bias_file
            header.add_history('Bias subtracted')
            prefix = 'b' + prefix

        return data, header, prefix

    @staticmethod
    def clean_column(_data, x0, y0, yf, n=5):
        t1 = _data[y0:yf, x0 - n:x0]
        t2 = _data[y0:yf, x0 + 1:x0 + n]
        t = np.hstack((t1, t2))
        _data[y0:yf, x0] = np.median(t)
        return _data

    def clean_columns(self, _data):
        bad_columns = [
            [167, 0, 512],
            [476, 0, 513],
            [602, 0, 513],
            [671, 0, 513],
            [810, 0, 513],
            [918, 0, 513],
            [917, 0, 513],
            [213, 513, 1024]
        ]

        for column in bad_columns:
            x0 = column[0]
            y0 = column[1]
            yf = column[2]
            _data = self.clean_column(_data, x0, y0, yf)
        return _data

    @staticmethod
    def clean_line(_data, x0, xf, y, n=5):
        t1 = _data[y - n:y, x0:xf]
        t2 = _data[y + 1:y + n, x0:xf]
        t = np.vstack((t1, t2))
        _data[y, x0:xf] = np.median(t)
        return _data

    def clean_lines(self, _data):
        bad_lines = [
            [214, 239, 688],
            [477, 516, 490],
            [387, 429, 455],
            [574, 603, 494],
            [574, 603, 493],
            [640, 672, 388],
            [604, 671, 388]
        ]

        for line in bad_lines:
            x0 = line[0]
            xf = line[1]
            y = line[2]
            _data = self.clean_line(_data, x0, xf, y)
        return _data

    @staticmethod
    def dark_subtraction(data, header, prefix, dark_file):

        if dark_file is not None:
            dark = pyfits.getdata(dark_file)
            data -= dark
            header['DARKFILE'] = dark_file
            prefix = 'd' + prefix
            header.add_history('Dark subtracted')

        return data, header, prefix

    @staticmethod
    def divide_by_flat(data, header, prefix, flat_file):

        if flat_file is not None:
            flat = pyfits.getdata(flat_file)
            data /= flat
            header['FLATFILE'] = flat_file
            header.add_history('Flat normalized')
            prefix = 'f' + prefix

        return data, header, prefix

    @staticmethod
    def divide_by_exposuretime(data, header, prefix, time):

        if time is True:
            try:
                exptime = float(header['EXPTIME'])
                data /= exptime
                header['UNITS'] = 'COUNTS/s'
                header.add_history('Divided by exposure time.')
                prefix = 't' + prefix
            except KeyError:
                pass

        return data, header, prefix

    @staticmethod
    def get_header(filename):

        fits_file = pyfits.open(filename)
        h0 = fits_file[0].header
        # noinspection PyUnusedLocal
        h1 = fits_file[1].header

        # TODO - Multiple header inheritance.
        # If there's any card that should be passed from the extentions
        # header to the main header, uncomment bellow.
        # for key in h1:
        #     if key not in h0:
        #         h0.set(key, value=h1[key], comment=h1.comments[key])

        h0.append('UNITS')
        h0.set('UNITS', value='COUNTS', comment='Pixel intensity units.')

        return h0

    @staticmethod
    def get_joined_data(filename):

        fits_file = pyfits.open(filename)
        w, h = str2pixels(fits_file[1].header['DETSIZE'])

        log.info(' > %s' % filename)

        # Correct for binning
        bin_size = np.array(fits_file[1].header['CCDSUM'].split(' '),
                            dtype=int)
        bw, bh = w[1] // bin_size[0], h[1] // bin_size[1]

        # Create empty full frame
        new_data = np.empty((bh, bw), dtype=float)

        # Process each extension
        for i in range(1, 5):
            tx, ty = str2pixels(fits_file[i].header['TRIMSEC'])
            bx, by = str2pixels(fits_file[i].header['BIASSEC'])

            data = fits_file[i].data
            trim = data[ty[0] - 1:ty[1], tx[0] - 1:tx[1]]
            bias = data[by[0] - 1:by[1], bx[0] - 1:bx[1]]

            # Collapse the bias columns to a single column.
            bias = np.median(bias, axis=1)

            # Fit and remove OVERSCAN
            x = np.arange(bias.size) + 1
            bias_fit_pars = np.polyfit(x, bias, 2)  # Last par = inf
            bias_fit = np.polyval(bias_fit_pars, x)
            bias_fit = bias_fit.reshape((bias_fit.size, 1))
            bias_fit = np.repeat(bias_fit, trim.shape[1], axis=1)

            trim = trim - bias_fit
            dx, dy = str2pixels(fits_file[i].header['DETSEC'])
            dx, dy = dx // bin_size[0], dy // bin_size[1]
            new_data[dy[0]:dy[1], dx[0]:dx[1]] = trim

        return new_data

    def main(self, list_of_files, bias_file=None, clean=False,
             cosmic_rays=False, dark_file=None, flat_file=None,
             glow_file=None, time=False):

        self.print_header()
        log.info('Processing data')
        list_of_files = sorted(list_of_files)

        for filename in list_of_files:
            prefix = "xj"

            # Get joined data
            data = self.get_joined_data(filename)

            # Build header
            header = self.get_header(filename)

            # BIAS subtraction
            data, header, prefix = self.bias_subtraction(
                data, header, prefix, bias_file
            )

            # DARK subtraction
            data, header, prefix = self.dark_subtraction(
                data, header, prefix, dark_file
            )

            # Remove cosmic rays and hot pixels
            data, header, prefix = self.remove_cosmic_rays(
                data, header, prefix, cosmic_rays
            )

            # Remove lateral glows
            data, header, prefix = self.remove_glows(
                data, header, prefix, glow_file
            )

            # FLAT division
            data, header, prefix = self.divide_by_flat(
                data, header, prefix, flat_file
            )

            # Normalize by the EXPOSURE TIME
            data, header, prefix = self.divide_by_exposuretime(
                data, header, prefix, time
            )

            # Removing bad column and line
            data = self.remove_central_bad_columns(data)

            # Clean known bad columns and lines
            data, header, prefix = self.clean_hot_columns_and_lines(
                data, header, prefix, clean
            )

            # Writing file
            header.add_history('Extensions joined using "sami_xjoin"')
            path, filename = os.path.split(filename)
            pyfits.writeto(os.path.join(path, prefix + filename), data,
                           header, clobber=True)

            log.info("\n All done!")

    @staticmethod
    def print_header():
        msg = (
            "\n SAMI - Join Extensions"
            " by Bruno Quint (bquint@astro.iag.usp.br)"
            " Mar 2015 - Version 0.4"
            "\n Starting program.")
        log.info(msg)

    @staticmethod
    def remove_cosmic_rays(data, header, prefix, cosmic_rays):

        if cosmic_rays:
            c = CosmicsImage(data.transpose(), gain=2.1, readnoise=1.83,
                             satlevel=65536)
            c.run(maxiter=4)

            data = c.cleanarray.transpose()

            header.add_history(
                'Cosmic rays and hot pixels removed using LACosmic')
            prefix += 'r'

        return data, header, prefix

    def remove_glows(self, data, header, prefix, glow_file):

        if glow_file is not None:
            # Create four different regions.
            regions = [
                [np.median(data[539:589, 6:56]),  # Top Left
                 np.median(data[539:589, 975:1019])],  # Top Right
                [np.median(data[449:506, 6:56]),  # Bottom Left
                 np.median(data[449:506, 975:1019])]  # Bottom Right
            ]
            min_std_region = np.argmin(regions) % 2
            # print(min_std_region)

            # The upper reg has background lower or equal to the lower reg
            midpt1 = regions[0][min_std_region]
            midpt2 = regions[1][min_std_region]
            diff = midpt2 - midpt1

            dark = pyfits.getdata(glow_file)
            dark = self.clean_columns(dark)
            dark = self.clean_lines(dark)

            dark_regions = [
                [np.median(dark[539:589, 6:56]),  # Top Left
                 np.median(dark[539:589, 975:1019])],  # Top Right
                [np.median(dark[449:506, 6:56]),  # Bottom Left
                 np.median(dark[449:506, 975:1019])]  # Bottom Right
            ]

            dark_midpt1 = dark_regions[0][min_std_region]
            dark_midpt2 = dark_regions[1][min_std_region]

            dark_diff = dark_midpt2 - dark_midpt1
            dark -= dark_midpt1

            k = diff / dark_diff
            temp_dark = dark * k
            data -= midpt1
            data -= temp_dark

            # print(k)

            header.add_history('Lateral glow removed using %s file' % glow_file)
            prefix = 'g' + prefix

        return data, header, prefix

    @staticmethod
    def set_debug(debug):

        if debug:
            log.basicConfig(level=log.DEBUG, format='%(message)s')

    @staticmethod
    def set_verbose(verbose):

        if verbose:
            log.basicConfig(level=log.INFO, format='%(message)s')
        else:
            log.basicConfig(level=log.WARNING, format='%(message)s')

    @staticmethod
    def remove_central_bad_columns(data):

        n_rows, n_columns = data.shape

        # Copy the central bad columns to a temp array
        temp_column = data[:, n_columns // 2 - 1:n_columns // 2 + 1]

        # Shift the whole image by two columns
        data[:, n_columns // 2 - 1:-2] = data[:, n_columns // 2 + 1:]

        # Copy the bad array in the end (right) of the image).
        data[:, -2:] = temp_column

        return data

    def clean_hot_columns_and_lines(self, data, header, prefix, clean):

        if clean is True:
            data = self.clean_columns(data)
            data = self.clean_lines(data)
            header.add_history('Cleaned bad columns and lines.')
            prefix = 'c' + prefix

        return data, header, prefix


# noinspection PyUnresolvedReferences
class CosmicsImage:
    """
    Docstring copied from the original cosmics.py file, obtained from
    U{http://obswww.unige.ch/~tewes/cosmics_dot_py/}

    About
    =====

    cosmics.py is a small and simple python module to detect and clean cosmic
    ray hits on images (numpy arrays or FITS), using scipy, and based on Pieter
    van Dokkum's L.A.Cosmic algorithm.

    L.A.Cosmic = Laplacian cosmic ray detection

    U{http://www.astro.yale.edu/dokkum/lacosmic/}

    (article : U{http://arxiv.org/abs/astro-ph/0108003})


    Additional features
    ===================

    I pimped this a bit to suit my needs :

    - Automatic recognition of saturated stars, including their full saturation
    trails. This avoids that such stars are treated as big cosmics. indeed
    saturated stars tend to get even uglier when you try to clean them. Plus
    they keep L.A.Cosmic iterations going on forever. This feature is mainly for
    pretty-image production. It is optional, requires one more parameter (a CCD
    saturation level in ADU), and uses some nicely robust morphology operations
    and object extraction.

    - Scipy image analysis allows to "label" the actual cosmic ray hits (i.e.
    group the pixels into local islands). A bit special, but I use this in the
    scope of visualizing a PSF construction.

    But otherwise the core is really a 1-to-1 implementation of L.A.Cosmic, and
    uses the same parameters. Only the conventions on how filters are applied at
     the image edges might be different.

    No surprise, this python module is much faster then the IRAF implementation,
    as it does not read/write every step to disk.

    Usage
    =====

    Everything is in the file cosmics.py, all you need to do is to import it.
    You need pyfits, numpy and scipy. See the demo scripts for example usages
    (the second demo uses f2n.py to make pngs, and thus also needs PIL).

    Your image should have clean borders, cut away prescan/overscan etc.


    Todo
    ====
    Ideas for future improvements :

    - Add something reliable to detect negative glitches (dust on CCD or small
      traps)
    - Top level functions to simply run all this on either numpy arrays or
      directly on FITS files
    - Reduce memory usage ... easy
    - Switch from signal to ndimage, homogenize mirror boundaries

    Malte Tewes, January 2010
    """

    def __init__(self, rawarray, pssl=0.0, gain=2.2, readnoise=10.0,
                 sigclip=5.0, sigfrac=0.3, objlim=5.0, satlevel=50000.0):
        """

        sigclip : increase this if you detect cosmics where there are none.
                  Default is 5.0, a good value for earth-bound images.

        objlim : increase this if normal stars are detected as cosmics. Default
                 is 5.0, a good value for earth-bound images.

        Constructor of the cosmic class, takes a 2D numpy array of your image as
        main argument.

        sigclip : laplacian-to-noise limit for cosmic ray detection

        objlim : minimum contrast between laplacian image and fine structure
                 image. Use 5.0 if your image is undersampled, HST, ...

        satlevel : if we find agglomerations of pixels above this level, we
                   consider it to be a saturated star and do not try to correct
                   and pixels around it. A negative satlevel skips this feature.

        pssl is the previously subtracted sky level !

        real   gain    = 1.8  # gain (electrons/ADU)	(0=unknown)
        real   readn   = 6.5  # read noise (electrons) (0=unknown)
        ##gain0  string statsec = "*,*" # section to use for automatic
                                        # computation of gain
        real   skyval  = 0.   # sky level that has been subtracted (ADU)
        real   sigclip = 3.0  # detection limit for cosmic rays (sigma)
        real   sigfrac = 0.5  # fractional detection limit for neighbouring
                              # pixels
        real   objlim  = 3.0  # contrast limit between CR and underlying object
        int    niter   = 1    # maximum number of iterations
        """

        # internally, we will always work "with sky".
        self.rawarray = rawarray + pssl

        # In lacosmiciteration() we work on this guy
        self.cleanarray = self.rawarray.copy()

        # All False, no cosmics yet
        self.mask = np.cast['bool'](np.zeros(self.rawarray.shape))

        self.gain = gain
        self.readnoise = readnoise
        self.sigclip = sigclip
        self.objlim = objlim
        self.sigcliplow = sigclip * sigfrac
        self.satlevel = satlevel

        self.pssl = pssl

        self.backgroundlevel = None  # only calculated and used if required.

        # a mask of the saturated stars, only calculated if required
        self.satstars = None

    def __str__(self):
        """
        Gives a summary of the current state, including the number of cosmic
        pixels in the mask etc.
        """
        stringlist = [
            "Input array : (%i, %i), %s" % (
                self.rawarray.shape[0], self.rawarray.shape[1],
                self.rawarray.dtype.name),
            "Current cosmic ray mask : %i pixels" % np.sum(self.mask)
        ]

        if self.pssl != 0.0:
            stringlist.append(
                "Using a previously subtracted sky level of %f" % self.pssl)

        if self.satstars is not None:
            stringlist.append(
                "Saturated star mask : %i pixels" % np.sum(self.satstars))

        return "\n".join(stringlist)

    def labelmask(self):
        """
        Finds and labels the cosmic "islands" and returns a list of dicts
        containing their positions. This is made on purpose for visualizations a
         la f2n.drawstarslist, but could be useful anyway.
        """

        log.debug("Labeling mask pixels ...")
        # We morphologicaly dilate the mask to generously connect "sparse"
        # cosmics :
        # dilstruct = np.ones((5,5))
        dilmask = ndimage.morphology.binary_dilation(self.mask,
                                                     structure=dilstruct,
                                                     iterations=1, mask=None,
                                                     output=None,
                                                     border_value=0, origin=0,
                                                     brute_force=False)
        # origin = 0 means center
        (labels, n) = ndimage.measurements.label(dilmask)

        # print "Number of cosmic ray hits : %i" % n
        # tofits(labels, "labels.fits", verbose = False)
        slicecouplelist = ndimage.measurements.find_objects(labels)

        # Now we have a huge list of couples of numpy slice objects giving a
        # frame around each object
        # For plotting purposes, we want to transform this into the center of
        # each object.
        if len(slicecouplelist) != n:
            # This never happened, but you never know ...
            raise RuntimeError("Mega error in labelmask !")
        centers = [[(tup[0].start + tup[0].stop) / 2.0,
                    (tup[1].start + tup[1].stop) / 2.0] for tup in
                   slicecouplelist]
        # We also want to know how many pixels where affected by each cosmic
        # ray. # Why ? Dunno... it's fun and available in scipy :-)
        sizes = ndimage.measurements.sum(self.mask.ravel(), labels.ravel(),
                                         np.arange(1, n + 1, 1))
        retdictlist = [{"name": "%i" % size, "x": center[0], "y": center[1]} for
                       (size, center) in zip(sizes, centers)]

        log.debug("Labeling done")

        return retdictlist

    def getdilatedmask(self, size=3):
        """
        Returns a morphologically dilated copy of the current mask.
        size = 3 or 5 decides how to dilate.
        """
        if size == 3:
            dilmask = ndimage.morphology.binary_dilation(self.mask,
                                                         structure=growkernel,
                                                         iterations=1,
                                                         mask=None, output=None,
                                                         border_value=0,
                                                         origin=0,
                                                         brute_force=False)
        elif size == 5:
            dilmask = ndimage.morphology.binary_dilation(self.mask,
                                                         structure=dilstruct,
                                                         iterations=1,
                                                         mask=None, output=None,
                                                         border_value=0,
                                                         origin=0,
                                                         brute_force=False)
        else:
            dilmask = self.mask.copy()

        return dilmask

    def clean(self, mask=None):
        """
        Given the mask, we replace the actual problematic pixels with the masked
        5x5 median value. This mimics what is done in L.A.Cosmic, but it's a
        bit harder to do in python, as there is no readymade masked median. So
        for now we do a loop...

        Saturated stars, if calculated, are also masked : they are not
        "cleaned", but their pixels are not used for the interpolation.

        We will directly change self.cleanimage. Instead of using the self.mask,
        you can supply your own mask as argument. This might be useful to apply
        this cleaning function iteratively. But for the true L.A.Cosmic, we
        don't use this, i.e. we use the full mask at each iteration.
        """

        if mask is None:
            mask = self.mask

        log.debug("Cleaning cosmic affected pixels ...")

        # So... mask is a 2D array containing False and True, where True means
        # "here is a cosmic".
        # We want to loop through these cosmics one by one.
        cosmicindices = np.argwhere(mask)
        # This is a list of the indices of cosmic affected pixels.
        # print cosmicindices

        # We put cosmic ray pixels to np.Inf to flag them :
        self.cleanarray[mask] = np.Inf

        # Now we want to have a 2 pixel frame of Inf padding around our image.
        w = self.cleanarray.shape[0]
        h = self.cleanarray.shape[1]
        padarray = np.zeros((w + 4, h + 4)) + np.Inf

        # that copy is important, we need 2 independent arrays
        padarray[2:w + 2, 2:h + 2] = self.cleanarray.copy()

        # The medians will be evaluated in this padarray, skipping the np.Inf.
        # Now in this copy called padarray, we also put the saturated stars to
        # np.Inf, if available :
        if self.satstars is not None:
            padarray[2:w + 2, 2:h + 2][self.satstars] = np.Inf
        # Viva python, I tested this one, it works...

        # A loop through every cosmic pixel :
        for cosmicpos in cosmicindices:

            x = cosmicpos[0]
            y = cosmicpos[1]

            # remember the shift due to the padding !
            cutout = padarray[x:x + 5, :y + 5].ravel()

            # print cutout
            # Now we have our 25 pixels, some of them are np.Inf, and we want to
            #  take the median
            goodcutout = cutout[cutout != np.Inf]
            # print np.alen(goodcutout)

            if np.alen(goodcutout) >= 25:
                # This never happened, but you never know ...
                raise RuntimeError("Mega error in clean !")
            elif np.alen(goodcutout) > 0:
                replacementvalue = np.median(goodcutout)
            else:
                # i.e. no good pixels : Shit, a huge cosmic, we will have to
                # improvise ...
                log.debug("OH NO, I HAVE A HUUUUUUUGE COSMIC !!!!!")
                replacementvalue = self.guessbackgroundlevel()

            # We update the cleanarray,
            # but measure the medians in the padarray, so to not mix things
            # up...
            self.cleanarray[x, y] = replacementvalue

        # That's it.
        log.debug("Cleaning done")

        # FYI, that's how the LACosmic cleaning looks in iraf :
        # noinspection PyPep8
        """
                imarith(outmask,"+",finalsel,outmask)
                imreplace(outmask,1,lower=1,upper=INDEF) # ok so outmask = 1 are the cosmics
                imcalc(outmask,inputmask,"(1.-10000.*im1)",verb-)
                imarith(oldoutput,"*",inputmask,inputmask)
                median(inputmask,med5,5,5,zloreject=-9999,zhi=INDEF,verb-)
                imarith(outmask,"*",med5,med5)
                if (i>1) imdel(output)
                imcalc(oldoutput//","//outmask//","//med5,output,"(1.-im2)*im1+im3",verb-)

                # =

                merging to full mask
                inputmask = 1.0 - 10000.0 * finalsel # So this is 1.0, but cosmics are very negative
                inputmask = oldoutput * inputmask # orig image, with very negative cosmics
                med5 = median of inputmask, but rejecting these negative cosmics
                # i dunno how to do this in python -> had to do the loop
                med5 = finalsel * med5 # we keep only the cosmics of this median
                # actual replacement :
                output = (1.0 - outmask)*oldoutput + med5 # ok
                """

    def findsatstars(self):
        """
        Uses the satlevel to find saturated stars (not cosmics !), and puts the
        result as a mask in self.satstars. This can then be used to avoid these
        regions in cosmic detection and cleaning procedures. Slow ...
        """
        log.debug("Detecting saturated stars ...")

        # DETECTION
        satpixels = self.rawarray > self.satlevel  # the candidate pixels

        # We build a smoothed version of the image to look for large stars and
        # their support :
        m5 = ndimage.filters.median_filter(self.rawarray, size=5, mode='mirror')
        # We look where this is above half the satlevel
        largestruct = m5 > (self.satlevel / 2.0)
        # The rough locations of saturated stars are now :
        satstarscenters = np.logical_and(largestruct, satpixels)

        log.debug("Building mask of saturated stars ...")

        # BUILDING THE MASK
        # The subtility is that we want to include all saturated pixels
        # connected to these saturated stars...
        # I haven't found a better solution then the double loop

        # We dilate the satpixels alone, to ensure connectivity in glitchy
        # regions and to add a safety margin around them.
        # dilstruct = np.array([[0,1,0], [1,1,1], [0,1,0]])

        dilsatpixels = ndimage.morphology.binary_dilation(satpixels,
                                                          structure=dilstruct,
                                                          iterations=2,
                                                          mask=None,
                                                          output=None,
                                                          border_value=0,
                                                          origin=0,
                                                          brute_force=False)
        # It turns out it's better to think large and do 2 iterations...

        # We label these :
        (dilsatlabels, nsat) = ndimage.measurements.label(dilsatpixels)
        # tofits(dilsatlabels, "test.fits")

        log.debug("We have %i saturated stars." % nsat)

        # The ouput, False for now :
        outmask = np.zeros(self.rawarray.shape)

        # we go through the islands of saturated pixels
        for i in range(1, nsat + 1):
            thisisland = dilsatlabels == i  # gives us a boolean array
            # Does this intersect with satstarscenters ?
            overlap = np.logical_and(thisisland, satstarscenters)
            # we add thisisland to the mask
            if np.sum(overlap) > 0:
                outmask = np.logical_or(outmask, thisisland)

        self.satstars = np.cast['bool'](outmask)

        log.debug("Mask of saturated stars done")

    def getsatstars(self):
        """
        Returns the mask of saturated stars after finding them if not yet done.
        Intended mainly for external use.
        """
        if not self.satlevel > 0:
            raise RuntimeError("Cannot determine satstars: "
                               "you gave satlevel <= 0 !")
        if self.satstars is None:
            self.findsatstars()
        return self.satstars

    def getmask(self):
        return self.mask

    def getrawarray(self):
        """
        For external use only, as it returns the rawarray minus pssl !
        """
        return self.rawarray - self.pssl

    def getcleanarray(self):
        """
        For external use only, as it returns the cleanarray minus pssl !
        """
        return self.cleanarray - self.pssl

    def guessbackgroundlevel(self):
        """
        Estimates the background level. This could be used to fill pixels in
        large cosmics.
        """
        if self.backgroundlevel is None:
            self.backgroundlevel = np.median(self.rawarray.ravel())
        return self.backgroundlevel

    def lacosmiciteration(self):
        """
        Performs one iteration of the L.A.Cosmic algorithm. It operates on
        self.cleanarray, and afterwards updates self.mask by adding the newly
        detected cosmics to the existing self.mask. Cleaning is not made
        automatically ! You have to call clean() after each iteration.
        This way you can run it several times in a row to to L.A.Cosmic
        "iterations". See function lacosmic, that mimics the full iterative
        L.A.Cosmic algorithm.

        Returns a dict containing
            - niter : the number of cosmic pixels detected in this iteration
            - nnew : among these, how many were not yet in the mask
            - itermask : the mask of pixels detected in this iteration
            - newmask : the pixels detected that were not yet in the mask

        If findsatstars() was called, we exclude these regions from the search.

        """

        log.debug("Convolving image with Laplacian kernel ...")

        # We subsample, convolve, clip negative values, and rebin to original
        # size
        subsam = subsample(self.cleanarray)
        conved = signal.convolve2d(subsam, laplkernel, mode="same",
                                   boundary="symm")
        cliped = conved.clip(min=0.0)
        # cliped = np.abs(conved) # unfortunately this does not work to find
        # holes as well ...
        lplus = rebin2x2(cliped)

        log.debug("Creating noise model ...")

        # We build a custom noise map, so to compare the laplacian to
        m5 = ndimage.filters.median_filter(self.cleanarray, size=5,
                                           mode='mirror')
        # We keep this m5, as I will use it later for the interpolation.
        m5clipped = m5.clip(min=0.00001)  # As we will take the sqrt
        noise = (1.0 / self.gain) * np.sqrt(
            self.gain * m5clipped + self.readnoise * self.readnoise)

        log.debug("Calculating Laplacian signal to noise ratio ...")

        # Laplacian signal to noise ratio :
        s = lplus / (2.0 * noise)  # the 2.0 is from the 2x2 subsampling
        # This s is called sigmap in the original lacosmic.cl

        # We remove the large structures (s prime) :
        sp = s - ndimage.filters.median_filter(s, size=5, mode='mirror')

        log.debug("Selecting candidate cosmic rays ...")

        # Candidate cosmic rays (this will include stars + HII regions)
        candidates = sp > self.sigclip
        nbcandidates = np.sum(candidates)

        log.debug("  %5i candidate pixels" % nbcandidates)

        # At this stage we use the saturated stars to mask the candidates, if
        # available :
        if self.satstars is not None:
            log.debug("Masking saturated stars ...")
            candidates = np.logical_and(np.logical_not(self.satstars),
                                        candidates)
            nbcandidates = np.sum(candidates)

            log.debug("  %5i candidate pixels not part of saturated stars" %
                      nbcandidates)

        log.debug("Building fine structure image ...")

        # We build the fine structure image :
        m3 = ndimage.filters.median_filter(self.cleanarray, size=3,
                                           mode='mirror')
        m37 = ndimage.filters.median_filter(m3, size=7, mode='mirror')
        f = m3 - m37
        # In the article that's it, but in lacosmic.cl f is divided by the
        # noise...
        # Ok I understand why, it depends on if you use sp/f or L+/f as
        # criterion. There are some differences between the article and the iraf
        # implementation. So I will stick to the iraf implementation.
        f /= noise
        # as we will divide by f. like in the iraf version.
        f = f.clip(min=0.01)

        log.debug("Removing suspected compact bright objects ...")

        # Now we have our better selection of cosmics :
        cosmics = np.logical_and(candidates, sp / f > self.objlim)
        # Note the sp/f and not lplus/f ... due to the f = f/noise above.

        nbcosmics = np.sum(cosmics)

        log.debug("  %5i remaining candidate pixels" % nbcosmics)

        # What follows is a special treatment for neighbors, with more relaxed
        # constains.

        log.debug("Finding neighboring pixels affected by cosmic rays ...")

        # We grow these cosmics a first time to determine the immediate
        # neighborhod  :
        growcosmics = np.cast['bool'](
            signal.convolve2d(np.cast['float32'](cosmics), growkernel,
                              mode="same", boundary="symm"))

        # From this grown set, we keep those that have sp > sigmalim
        # so obviously not requiring sp/f > objlim, otherwise it would be
        # pointless
        growcosmics = np.logical_and(sp > self.sigclip, growcosmics)

        # Now we repeat this procedure, but lower the detection limit to
        # sigmalimlow :

        finalsel = np.cast['bool'](
            signal.convolve2d(np.cast['float32'](growcosmics), growkernel,
                              mode="same", boundary="symm"))
        finalsel = np.logical_and(sp > self.sigcliplow, finalsel)

        # Again, we have to kick out pixels on saturated stars :
        if self.satstars is not None:
            log.debug("Masking saturated stars ...")
            finalsel = np.logical_and(np.logical_not(self.satstars), finalsel)

        nbfinal = np.sum(finalsel)

        log.debug("  %5i pixels detected as cosmics" % nbfinal)

        # Now the replacement of the cosmics...
        # we outsource this to the function clean(), as for some purposes the
        # cleaning might not even be needed. Easy way without masking would be:
        # self.cleanarray[finalsel] = m5[finalsel]

        # We find how many cosmics are not yet known :
        newmask = np.logical_and(np.logical_not(self.mask), finalsel)
        nbnew = np.sum(newmask)

        # We update the mask with the cosmics we have found :
        self.mask = np.logical_or(self.mask, finalsel)

        # We return
        # (used by function lacosmic)

        return {"niter": nbfinal, "nnew": nbnew, "itermask": finalsel,
                "newmask": newmask}

    # noinspection PyUnusedLocal
    @staticmethod
    def findholes():
        """
        Detects "negative cosmics" in the cleanarray and adds them to the mask.
        This is not working yet.
        """
        # TODO Implement this
        pass

        # noinspection PyPep8
        """
                if verbose == None:
                    verbose = self.verbose

                if verbose :
                    print "Finding holes ..."

                m3 = ndimage.filters.median_filter(self.cleanarray, size=3, mode='mirror')
                h = (m3 - self.cleanarray).clip(min=0.0)

                tofits("h.fits", h)
                sys.exit()

                # The holes are the peaks in this image that are not stars

                #holes = h > 300
                """
        # noinspection PyPep8
        """
                subsam = subsample(self.cleanarray)
                conved = -signal.convolve2d(subsam, laplkernel, mode="same", boundary="symm")
                cliped = conved.clip(min=0.0)
                lplus = rebin2x2(conved)

                tofits("lplus.fits", lplus)

                 m5 = ndimage.filters.median_filter(self.cleanarray, size=5, mode='mirror')
                 m5clipped = m5.clip(min=0.00001)
                 noise = (1.0/self.gain) * np.sqrt(self.gain*m5clipped + self.readnoise*self.readnoise)

                 s = lplus / (2.0 * noise) # the 2.0 is from the 2x2 subsampling
                 # This s is called sigmap in the original lacosmic.cl

                 # We remove the large structures (s prime) :
                 sp = s - ndimage.filters.median_filter(s, size=5, mode='mirror')

                 holes = sp > self.sigclip
                """
        """
        # We have to kick out pixels on saturated stars :
        if self.satstars != None:
             if verbose:
                 print "Masking saturated stars ..."
             holes = np.logical_and(np.logical_not(self.satstars), holes)

        if verbose:
            print "%i hole pixels found" % np.sum(holes)

        # We update the mask with the holes we have found :
        self.mask = np.logical_or(self.mask, holes)
        """

    def run(self, maxiter=4):
        """
        Full artillery :-)
            - Find saturated stars
            - Run maxiter L.A.Cosmic iterations (stops if no more cosmics are
              found)

        Stops if no cosmics are found or if maxiter is reached.
        """

        if self.satlevel > 0 and self.satstars is None:
            self.findsatstars()

        log.debug("Starting %i L.A.Cosmic iterations ..." % maxiter)
        for i in range(1, maxiter + 1):
            log.debug("Iteration %i" % i)

            iterres = self.lacosmiciteration()
            log.debug("%i cosmic pixels (%i new)" %
                      (iterres["niter"], iterres["nnew"]))

            # self.clean(mask = iterres["mask"]) # No, we want clean to operate
            # on really clean pixels only! Thus we always apply it on the full
            # mask, as lacosmic does :
            self.clean()
            # But note that for huge cosmics, one might want to revise this.
            # Thats why I added a feature to skip saturated stars !

            if iterres["niter"] == 0:
                break


# noinspection PyUnusedLocal
def rebin(a, newshape):
    """
    Auxiliary function to rebin an ndarray a.
    U{http://www.scipy.org/Cookbook/Rebinning}

    a=rand(6,4); b=rebin(a,(3,2))
    """

    shape = a.shape
    len_shape = len(shape)
    factor = np.asarray(shape) / np.asarray(newshape)
    # print factor
    # noinspection PyPep8,PyPep8,PyPep8
    ev_list = ['a.reshape('] + \
              ['newshape[%d],factor[%d],' % (i, i) for i in xrange(len_shape)] + \
              [')'] + ['.sum(%d)' % (i + 1) for i in xrange(len_shape)] + \
              ['/factor[%d]' % i for i in xrange(len_shape)]

    return eval(''.join(ev_list))


def rebin2x2(a):
    """
    Wrapper around rebin that actually rebins 2 by 2.
    """
    inshape = np.array(a.shape)
    if not (inshape % 2 == np.zeros(
            2)).all():  # Modulo check to see if size is even
        raise RuntimeError("I want even image shapes !")

    return rebin(a, inshape / 2)


def subsample(a):
    """
    Piece of code extracted from cosmics.py.
    Original source obtained from:
    U{http://www.astro.yale.edu/dokkum/lacosmic/}

    Returns a 2x2-subsampled version of array a (no interpolation, just cutting
    pixels in 4).

    The version below is directly from the scipy cookbook on rebinning :
    U{http://www.scipy.org/Cookbook/Rebinning}

    There is
    ndimage.zoom(cutout.array, 2, order=0, prefilter=False),
    but it makes funny borders.
    """

    newshape = (2 * a.shape[0], 2 * a.shape[1])

    slices = [slice(0, old, float(old) / new) for old, new in
              zip(a.shape, newshape)]

    coordinates = np.mgrid[slices]

    # choose the biggest smaller integer index
    indices = coordinates.astype('i')

    return a[tuple(indices)]


def str2pixels(my_string):
    my_string = my_string.replace('[', '')
    my_string = my_string.replace(']', '')
    x, y = my_string.split(',')

    x = x.split(':')
    y = y.split(':')

    # "-1" fix from IDL to Python
    x = np.array(x, dtype=int)
    y = np.array(y, dtype=int)

    return x, y


if __name__ == '__main__':
    # Parsing Arguments ---
    parser = argparse.ArgumentParser(
        description="Join extensions existent in a single FITS file."
    )

    parser.add_argument('-b', '--bias', type=str, default=None,
                        help="Consider BIAS file for subtraction.")
    # TODO Enable Clean Hot Columns and Lines
    # parser.add_argument('-c', '--clean', action='store_true',
    #                     help="Clean known bad columns and lines by taking the "
    #                          "median value of their neighbours.")
    parser.add_argument('-d', '--dark', type=str, default=None,
                        help="Consider DARK file for subtraction.")
    parser.add_argument('-D', '--debug', action='store_true',
                        help="Turn on DEBUG mode (overwrite quiet mode).")
    parser.add_argument('-f', '--flat', type=str, default=None,
                        help="Consider FLAT file for division.")
    parser.add_argument('-g', '--glow', type=str, default=None,
                        help="Consider DARK file to correct lateral glows.")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run quietly.")
    # TODO Enable Cosmic Ray Removals
    # parser.add_argument('-r', '--rays', action='store_true',
    #                     help='Use LACosmic.py to remove cosmic rays and hot '
    #                          'pixels.')
    parser.add_argument('-t', '--exptime', action='store_true',
                        help="Divide by exposure time.")
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                        help="input filenames.")

    pargs = parser.parse_args()
    # SAMI_XJoin(
    #     pargs.files, bias_file=pargs.bias, clean=pargs.clean,
    #     cosmic_rays=pargs.rays, dark_file=pargs.dark, debug=pargs.debug,
    #     flat_file=pargs.flat, glow_file=pargs.glow, time=pargs.exptime,
    #     verbose=not pargs.quiet
    # )

    SAMI_XJoin(pargs.files, bias_file=pargs.bias, dark_file=pargs.dark,
        debug=pargs.debug, flat_file=pargs.flat, glow_file=pargs.glow,
        time=pargs.exptime, verbose=not pargs.quiet)
