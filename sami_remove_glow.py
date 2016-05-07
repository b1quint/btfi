#!/usr/bin/env python3
# -*- coding: utf8 -*-
import numpy as np
import pyfits
from scipy import ndimage

from glob import glob

print(np.__version__)

# dark_file = '/data/20160505/007/DARKscaled.fits'
#dark_file = '/data/20160505/007/DARKs120.fits'
dark_file = '/data/20160505/007/DARK300.fits'
#dark_file = '/data/20160506/008/DARK60s.fits'

data = pyfits.getdata(dark_file)
hdr = pyfits.getheader(dark_file)

#data = ndimage.median_filter(data, [1,5])
#data = ndimage.median_filter(data, [5,1])
#hdr.add_history("3x3 median filter applied")

def clean_columns(data, x0, y0, yf, n=5):
    t1 = data[y0:yf, x0-n:x0]
    t2 = data[y0:yf, x0+1:x0+n]
    t = np.hstack((t1, t2))
    print(y0, yf, x0-n, x0, x0+1, x0+n)
    data[y0:yf, x0] = np.median(t)
    return data

def clean_line(data, x0, xf, y, n=5):
    t1 = data[y-n:y, x0:xf]
    t2 = data[y+1:y+1, x0:xf]
    t = np.vstack((t1, t2))
    data[y, x0:xf] = np.median(t)
    return data

data = clean_columns(data, 167, 0, 512)
data = clean_columns(data, 476, 0, 513)
data = clean_columns(data, 602, 0, 513) 
data = clean_columns(data, 671, 0, 513)
data = clean_columns(data, 810, 0, 513)
data = clean_columns(data, 918, 0, 513)
data = clean_columns(data, 917, 0, 513)
data = clean_columns(data, 213, 513, 1024)

lines = [[214, 239, 688],
         [477, 516, 490],
         [387, 429, 455],
         [574, 603, 494],
         [574, 603, 493],
         [640, 672, 388],
         [604, 671, 388]
         ]
for line in lines:
    data = clean_line(data, line[0], line[1], line[2])


dark_midpt1 = np.median(data[539:589, 999:1009])
dark_midpt2 = np.median(data[449:506, 975:1019])
dark_diff = dark_midpt2 - dark_midpt1
data -= dark_midpt1
#pyfits.writeto(dark_file.replace('.fits', '.median.fits'), data, hdr)
print('Dark Stats')
print(dark_diff)

science_file = '/data/20160505/007/tbxj*.fits' 
files = glob(science_file)

for f in files:

    d = pyfits.getdata(f)
    h = pyfits.getheader(f)

    midpt1 = np.median(d[539:589, 999:1009])
    midpt2 = np.median(d[449:506, 975:1019])
    diff = midpt2 - midpt1
    print(diff)

    k = diff / dark_diff
    temp_dark = data * k
    d -= midpt1
    d -= temp_dark
    print("k = ", k)

    h.add_history('Lateral glow subtracted')
    pyfits.writeto(f.replace('/tbxj', '/ctbxj'), d, h, clobber=True)



