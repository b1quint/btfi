#/usr/bin/python
# -*- coding: utf8 -*-
"""
    ARtificial Imunologic ALgorithm applied to Spectra

    This script is a first draft of the algorithm described above. More
    information will be given during the code.

    By Bruno C. Quint (bquint@astro.iag.usp.br)
    Dec 2015
"""
from __future__ import division, print_function

import argparse
import astropy.io.fits as pyfits
import logging
import matplotlib.pyplot as plt
import numpy as np

from scipy import stats

__author__ = 'Bruno Quint'
__ver__ = '0.0.0'
__date__ = '2015.12.14'

log_formatter = logging.Formatter(" %(levelname)s - %(message)s")

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_formatter)

log = logging.getLogger()
log.addHandler(stream_handler)
log.setLevel(logging.DEBUG)

class Main(object):

    def __init__(self, verbose=False, debug=False):
        """
        Main class constructor.

        @param verbose: Set verbose level on.
        @param debug: Set debug level on. This overwrites verbose parameter.
        @return None
        """
        # Setting verbose/debug levels
        if verbose is False:
            log.setLevel(logging.ERROR)
        elif verbose is True:
            log.setLevel(logging.INFO)

        if debug is True:
            log.setLevel(logging.DEBUG)

        return

    def gaussian(self, x, x0, fwhm):
        """
        @param x: 1D array containing X values.
        @param x0: center peak position.
        @param fwhm: gaussian's fwhm.
        @return y: gaussian shape with peak equals to 1.
        """

        A = 1.0
        B = x0
        C = fwhm / (2 * np.sqrt(2 * np.log(2)))

        y = A * np.exp(-(x - B) ** 2 / (2 * C ** 2))

        return y

    def sim_single_peak(self, size=50, fwhm=3, snr=10.):
        """
        Simulates data for debugging. It basically creates a fake spectrum with
         a given SNR.

        @param size: array's size.
        @param fwhm: peak's FWHM.
        @param snr: signal-to-noise level.
        @return x, y: two arrays containing X and Y values for gaussian.
        """

        # Create empty data.
        data = np.zeros(size)

        # Get random peak position
        peak_pos = np.random.random_integers(0, size)
        x = np.arange(size)

        # Create gaussian signal
        y = self.gaussian(x, peak_pos, fwhm)

        # Add random noise
        n = np.random.rand(size)
        y = snr * y + n

        return x, y

    def run(self, filename=None, pop=5, gen=100, wsize=5):
        """
        Run main thread.

        @param filename: Input FITS 1D spectrum.
        @param pop: number of antibodies
        @param gen: number of generations (interactions)
        @param wsize: window size.
        @return None.
        """

        # Print a cool header for the user
        log.info("ARIALS - ARtificial Imunologic ALgorithm applied to Spectra")
        log.info("by Bruno C. Quint - bquint at ctio.noao.edu")
        log.info("Version %s" % __ver__)
        log.info("Date %s" % __date__)

        # Loading file/Simulating data
        if filename is None:
            x, y = self.sim_single_peak(snr=5)
        else:
            y = pyfits.getdata(filename)
            h = pyfits.getheader(filename)
            x = np.arange(y.size)

            try:
                x = (x - h['CRPIX1'] + 1) * h['CDELT1'] + h['CRVAL1']
            except KeyError:
                pass

        # Start ARIALS
        peak = Arials(x, y, pop, gen, wsize).run()
        log.info("Peak found at %02d" % peak)

        return


class Arials(object):

    def __init__(self, x_data, y_data, population, generation, sub_region_size):
        """
        Class constructor that will receive X and Y data.
        @param x_data: numpy.1darray with X values.
        @param y_data: numpy.1darray with Y values.
        @param population: number of anti-bodies.
        @param generation: number of generations.
        @param sub_region_size: number of sub-regions.
        @return None
        """
        self.x = x_data
        self.y = y_data
        self.population = population
        self.generation = generation
        self.sub_region_size = sub_region_size
        return

    @staticmethod
    def get_antibody(size, w_size):
        """
        Return the position of a single antibody.

        @param size: an intteger containing the spectrum size.
        @param w_size: the size of the window.
        """
        x = np.random.random_integers(1, size / w_size - 1)
        return x

    def get_antibodies(self, size, w_size, n_bodies):
        """
        Returns a list containing a set of antibodies.

        @param size: an intteger containing the spectrum size.
        @param wSize: the size of the windows
        @param nBodies: the number of antibodies
        """
        abody = []
        for i in range(n_bodies):
            abody.append(self.get_antibody(size, w_size))

        return abody

    @staticmethod
    def get_clone(antibody):
        """
        Returns a random antibody clone.

        @param antibody: is a list containing the coordinates of a single
        antibody.
        """
        v = range(9)
        v.remove(4)
        i = np.random.random_integers(0, 7)
        v = v[i]
        x = (v / 3) - 1
        y = (v % 3) - 1
        return np.array([antibody + x])

    @staticmethod
    def get_region_sum(y, i, w):
        """
        Returns the sum of the values of the pixels in a windows with 'w' size
        around the channel i.

        @param y: data being analysed.
        @param i: center of the region to be summed.
        @param w: window size.
        """
        x0 = max(i - w // 2, 0)
        xf = min(i + w // 2, y.size - 1)

        if x0 < 0:
            log.warning(" Warning: window out of image (line < 0)")

        if xf > y.shape[0]:
            log.warning(" Warning: window out of image (line > image size)")

        return y[x0:xf].sum()

    def run(self, show=True):
        """
        Run ARIALS algorithm.
        @param show: show plots.
        @return None.
        """
        x = self.x
        y = self.y
        gen = self.generation
        pop = self.population
        srs = self.sub_region_size

        # Show data
        if show:
            fig = plt.figure()
            ax = plt.axes(xlim=(x.min(), x.max()), ylim=(y.min(), y.max()))
            plt.plot(x, y, 'k-', lw=2)
            plt.grid()

        # Generate antibodies.
        antibodies = self.get_antibodies(x.size, 1, pop)

        # First generation
        my_sum = np.zeros(len(antibodies))

        # Evolving
        for n in range(gen):
            for i in range(pop):

                # Interact over an antibody
                ab = antibodies[i]
                my_sum[i] = self.get_region_sum(y, ab, srs)

                # Create clones
                clone1 = max(ab - 2, 0)
                clone2 = min(ab + 2, y.size - 1)

                # Sum regions near clones
                sum1 = self.get_region_sum(y, clone1, srs)
                sum2 = self.get_region_sum(y, clone2, srs)

                # Analyse regions
                sumf = sum1 if sum1 > sum2 else sum2
                clonef = clone1 if sum1 > sum2 else clone2
                ab = ab if my_sum[i] > sumf else clonef

                log.debug("GEN %03d - ANB %02d at %02d" % (n, i, ab))

                if show:
                    plt.plot(ab, y[ab], alpha=1, mec='None', mfc='red', marker='o')

            new_antibody = self.get_antibody(y.size, srs)
            for j in range(self.population // 2):
                antibodies[my_sum.argsort()[j]] = new_antibody

        # Show data
        if show:
            s = ''
            for ab in antibodies:
                plt.axvline(ab, c='r')
                s += '\t %02d' % ab
            log.info(s)
            plt.show()

        peak = stats.mode(antibodies)[0].astype(int)
        return x[peak]


if __name__ == '__main__':

    # Parse arguments ---------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="ARIALS - ARtificial Immunologic ALgorithm applied to " +
                    "Spectra")

    parser.add_argument('-f', '--filename', type=str, default=None,
                        help="Input spectrum.")

    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Turn verbose on.")

    parser.add_argument('-D', '--DEBUG', action='store_true',
                        help="Turn DEBUG on.")

    parser.add_argument('-p', '--population', type=int, default=10,
                        help='Number of antibodies present.')

    parser.add_argument('-g', '--generation', type=int, default=100,
                        help='Number of generations.')

    parser.add_argument('-w', '--windowsize', type=int, default=3,
                        help='Window size.')

    args = parser.parse_args()

    main = Main(verbose=True, debug=args.DEBUG)
    main.run(args.filename, pop=args.population,
             gen=args.generation, wsize=args.windowsize)
