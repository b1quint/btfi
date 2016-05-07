#!/usr/bin/env python3
# -*- coding: utf8 -*-
import numpy as np
import pyfits
from scipy import ndimage

from glob import glob


columns =[
        [167, 0, 512], 
        [476, 0, 513], 
        [602, 0, 513],   
        [671, 0, 513], 
        [810, 0, 513], 
        [918, 0, 513], 
        [917, 0, 513], 
        [213, 513, 1024]
        ]

lines = [[214, 239, 688],
         [477, 516, 490],
         [387, 429, 455],
         [574, 603, 494],
         [574, 603, 493],
         [640, 672, 388],
         [604, 671, 388]
         ]        
        
def main():

    filename = '/data/20160505/007/HD303056.fits' 
    global columns
    global lines
    
    data = pyfits.getdata(filename)
    hdr = pyfits.getheader(filename)
    
    for col in columns:
        data = clean_columns(data, col[0], col[1], col[2])
    
    for line in lines:
        data = clean_line(data, line[0], line[1], line[2])
       
    pyfits.writeto(filename.replace('.fits', '.clean.fits'), 
                    data, hdr, clobber=True)
    

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

if __name__ == '__main__':
    main()
