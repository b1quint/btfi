#/usr/bin/python
# -*- coding: utf8 -*-
"""
    SAMI - Remove Foorprint

    This script uses dark files to analyse model and subtract lateral footprints
    present on SAMI's images.
"""
from __future__ import division, print_function

import argparse
import astropy.io.fits as pyfits
import logging
import matplotlib.pyplot as plt

from astropy.visualization import mpl_normalize, stretch

__author__ = 'Bruno Quint'
__ver__ = '0.0.0'
__date__ = '2015.12.14'

log_formatter = logging.Formatter(" %(levelname)s - %(message)s")

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_formatter)

log = logging.getLogger()
log.addHandler(stream_handler)
log.setLevel(logging.ERROR)


class Main(object):

    def __init__(self, dark_file, list_of_files, verbose_level=2):

        try:
            assert verbose_level in [0, 1, 2]
        except AssertionError:
            import sys
            log.error("Invalid verbose level. Leaving now.")
            sys.exit(1)

        self.dark_file = dark_file
        self.list_of_files = list_of_files

        if verbose_level == 0:
            log.setLevel(logging.ERROR)
        elif verbose_level == 1:
            log.setLevel(logging.INFO)
        else:
            log.setLevel(logging.DEBUG)

        return

    def run(self):

        # Print a cool header for the user ---
        log.info("SAMI Remove Footprint")
        log.info("by Bruno C. Quint - bquint at ctio.noao.edu")
        log.info("Version %s" % __ver__)
        log.info("Date %s" % __date__)

        # Open the dark file ---
        dark_data = pyfits.getdata(self.dark_file)




        return

    def show_areas(self, data):
        """
        Show where all the calculations are being executed.
        """
        config = {'interpolation': 'nearest',
                  'origin': 'lower',
                  'cmap': 'hot'}

        vmin = data.mean() - 0.0 * data.std()
        vmax = data.mean() + 0.1 * data.std()

        plt.imshow(data, vmin=vmin, vmax=vmax, **config)
        plt.grid(c='white')

        ax = plt.gca()
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])

        plt.show()

        return



class bcolors:
    """
    Colors to be used for a better interaction with the user. Source from:
    http://stackoverflow.com/a/287944/2333908
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


if __name__ == '__main__':

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(
        description="Uses a master dark file to model and subtract the lateral"+
                    " footprints sometimes present on SAMI's images."
    )

    parser.add_argument('dark_file', type=str, help="Dark file.")
    parser.add_argument('files', type=str, nargs='+', help="input filenames.")
    parser.add_argument('-v','--verbose', type=int, default=0,
                        help="Verbose level: 0 - None, 1 - Info, 2 - Debug.")

    args = parser.parse_args()

    # Main thread ---
    main = Main(args.dark_file, args.files, verbose_level=args.verbose)
    main.run()
