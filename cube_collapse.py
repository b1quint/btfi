#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    Cube Collapse

    As its title says, this file collapses a data-cube in the spectral (Z)
    direction using by summing each frame, using the mean or the median.
    The default operation used is the average.

    By Bruno C. Quint - bquint@ctio.noao.edu
    2015.10.23
"""
from __future__ import division, print_function

import argparse
import astropy.io.fits as pyfits
import logging

__author__ = 'Bruno Quint'
__version__= '0.0.0'


class Main(object):

    def __init__(self):
        return

    def parse_arguments(self):
        """Parse the arguments given by the user on the commmand line."""
        parser = argparse.ArgumentParser(description="Collapse a data-cube.")
        parser.add_argument('input_filename', type=str, nargs=1,
                            help="Input filename.")

        args = parser.parse_args()
        return args

    def run(self):
        return

if __name__ == '__main__':
    main = Main()
    main.run()
