#!/usr/bin/env python2
import sys
sys.path.append("../.") # The directory that contains cosmic.py

import astropy.io.fits as pyfits
import cosmics
import numpy as np
import os

# File information
path = '/home/bquint/Data/20151203/012'
# name = 'NGC2440_3D_014_phc_zcenter.fits'
name = 'NGC2440_3D_015_phc_zcenter.fits'
name = 'xjfp_sami_C001.602.fits'
# name = 'temp.fits'

filename = os.path.join(path, name)

# Read the FITS :
cube = pyfits.getdata(filename)
header = pyfits.getheader(filename)

cube = np.asarray(cube)
print(cube.shape)

# array is a 3D numpy array
for i in range(cube.shape[0]):

    channel = cube[i].transpose()
    print(i, channel.shape)

    # Build the object :
    c = cosmics.cosmicsimage(channel, gain=2.1, readnoise=10.0, sigclip=3.0,
                             sigfrac=0.3, objlim=5.0)

    # Run the full artillery :
    c.run(maxiter=4)

    clean_channel = c.cleanarray.transpose()
    print(clean_channel.shape)
    cube[i] = clean_channel

print(cube.shape)
pyfits.writeto(filename.replace('.fits', '_crclean.fits'), cube, header, clobber=True)
