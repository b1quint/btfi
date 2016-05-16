# The following two lines are only needed as cosmic.py is not in this directory nor in the python path.
# They would not be required if you copy cosmics.py in this directory.
import sys
import os
import pyfits
sys.path.append("../.") # The directory that contains cosmic.py

import cosmics

# File information
path = '/home/bquint/Data/20151203/012'
name = 'xjfp_sami_C001.602.fits'
filename = os.path.join(path, name)

# Read the FITS :
array = pyfits.getdata(filename)

# Build the object :
c = cosmics.cosmicsimage(array, gain=2.1, readnoise=10.0, sigclip = 3.0, sigfrac = 0.3, objlim = 5.0)
# There are other options, check the manual...

# Run the full artillery :
c.run(maxiter = 4)

array = c.cleanarray
pyfits.writeto(filename.replace('.fits', '_crclean.fits'), array, clobber=True)
