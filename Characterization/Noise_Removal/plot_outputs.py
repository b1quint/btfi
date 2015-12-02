#/usr/bin/python
# -*- coding: utf8 -*-
"""
    Standard Deviation Vs the number of BIAS combined - Plot

    by Bruno C. Quint
"""
from __future__ import division, print_function

import glob
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'Bruno Quint'

def main():

    # Two subplots, the axes array is 1-d
    f, axarr = plt.subplots(2, sharex=True)

    files = glob.glob('*.dat')
    for f in files:
        n, a, s = np.loadtxt(f, unpack=True)
        axarr[0].plot(n, a, label=f)
        axarr[0].set_title('N x Average')
        axarr[0].grid()

        axarr[1].plot(n, s, label=f)
        axarr[1].set_title('N x STD')
        axarr[1].legend(loc='best')

    plt.show()
    print('.')


if __name__ == '__main__':
    main()
