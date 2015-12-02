#/usr/bin/python
# -*- coding: utf8 -*-
"""
    Standar Deviation Vs the number of BIAS combined

    by Bruno C. Quint
"""
from __future__ import division, print_function

import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pyraf.iraf as iraf
import glob
import os
import random

__author__ = 'Bruno Quint'

def main():

    config = {'combine': 'average',
              'reject': 'minmax',
              'nlow': 3,
              'nhigh': 3}

    path = '/btfidr/data/20130902/007/cam2'
    kernel = 'biasEMOFF_A*'
    list_of_files = glob.glob(os.path.join(path, kernel))
    random.shuffle(list_of_files)
    number_of_files = len(list_of_files)

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)

    n, a, s = [], [], []
    for l in range(2, number_of_files):

        temp_list = 'list'
        buffer = open(temp_list, mode='w')
        for i in range(l):
            buffer.write(list_of_files[i] + "\n")
        buffer.close()

        iraf.imcombine(input='@' + temp_list,
                       output='temp_output.fits',
                       **config)

        data = pyfits.getdata('temp_output.fits')[700:899,700:899]

        n.append(l)
        a.append(np.average(data))
        s.append(np.std(data))

        os.remove('temp_output.fits')
        os.remove(temp_list)

    output = 'output.dat'
    buffer = open(output, mode='w')

    for key in config.keys():
        buffer.write('# {}: {}\n'.format(key, config[key]))

    buffer.write('# N\tAVG\tSTD\n')
    for i in range(len(n)):
        buffer.write(u'{0:3d}\t{1:5.1f}\t{2:5.1f}\n'.format(n[i], a[i], s[i]))

    buffer.close()
    return 0


if __name__ == '__main__':
    main()
