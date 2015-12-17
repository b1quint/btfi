# IPython log file

get_ipython().magic('logstart')
import astropy.io.fits as pyfits
data = pyfits.getdata('bquint_datared/cal6600-20_013--phc-residual.fits')
hdr = pyfits.getheader('bquint_datared/cal6600-20_013--phc-residual.fits')
import np
import numpy as np
z = np.arange(data.shape[0])
z
center = (data * z).sum(axis=1).axis(2) / data.sum()
center = (data * z).sum(axis=1).axis(2) / data.sum()
z = z.reshape((z.size, 1, 1))
z = z.repeat(1028, axis=1)
z = z.repeat(1024, axis=2)
z.shape
data.shape
center  = (data * z).sum(axis=2).sum(axis=1) / data.sum(axis=2).sum(axis=1)
z = (z - hdr['CRPIX3'] + 1) * hdr['C3_3'] + hdr['CRVAL3']
center  = (data * z).sum(axis=2).sum(axis=1) / data.sum(axis=2).sum(axis=1)
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--center.fits')
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--center.fits', data, header)
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--center.fits', data, hdr)
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--center.fits', center, hdr, clobber=True)
center.shape
center  = (data * z).sum(axis=0) / data.sum(axis=0)
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--center.fits', center, hdr, clobber=True)
import matplotlib.pyplot as plt
plt.hist(center.ravel())
plt.show()
plt.ion()
plt.hist(center.ravel(), bins=100)
plt.savefig('bquint_datared/cal6600-20_013--phc-residual--center--hist.png')
fwhm = np.sqrt((data * z * z).sum(axis=0) / data.sum(axis=0) - center ** 2)
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--fwhm.fits', fwhm, hdr, clobber=True)
fwhm = np.sqrt((data * z * z).sum(axis=0) / data.sum(axis=0) - center ** 2) / (2 * np.sqrt(2 * np.ln(2))
)
fwhm = np.sqrt((data * z * z).sum(axis=0) / data.sum(axis=0) - center ** 2) / (2 * np.sqrt(2 * np.log(2)))
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--fwhm.fits', fwhm, hdr, clobber=True)
https://en.wikipedia.org/wiki/Image_moment
def lorentz(p,x):
        return p[2] / (1.0 + (x / p[0] - 1.0)**4 * p[1]**2)

def lorentz(p,x):
        return p[2] / (1.0 + (x / p[0] - 1.0)**4 * p[1]**2)

from scipy.optimize import leastsq
def lorentz(p,x):
        return (p[0] / (np.pi * p[2])) * (1 / (1 + ((x - p[1])/p[2])**2))

def errorfunc(p,x,zz):
    return lorentz(p,x)-zz

for i in range(data.shape[2]):
    for j in range(data.shape[1]):
        p[0] = data[:,j,i].max()
        p[1] = np.argmax(data[:,j,i])
        p[2] = fwhm[j,i]
        solp, ier = leastsq(errorfunc, p)
        fwhm[j,i] = solp[2]
        center[j,i] = solp[1]
        
p = [0, 0, 0]
for i in range(data.shape[2]):
    for j in range(data.shape[1]):
        p[0] = data[:,j,i].max()
        p[1] = np.argmax(data[:,j,i])
        p[2] = fwhm[j,i]
        solp, ier = leastsq(errorfunc, p)
        fwhm[j,i] = solp[2]
        center[j,i] = solp[1]
        
for i in range(data.shape[2]):
    for j in range(data.shape[1]):
        p[0] = data[:,j,i].max()
        p[1] = np.argmax(data[:,j,i])
        p[2] = fwhm[j,i]
        solp, ier = leastsq(errorfunc, p, args=(z[:,j,i],data[:,j,i]))
        fwhm[j,i] = solp[2]
        center[j,i] = solp[1]
        print(i,j,p,solp)
        
for i in range(data.shape[2]):
    for j in range(data.shape[1]):
        p[0] = data[:,j,i].max()
        p[1] = np.argmax(data[:,j,i])
        p[2] = fwhm[j,i]
        solp, ier = leastsq(errorfunc, p, args=(z[:,j,i],data[:,j,i]))
        fwhm[j,i] = solp[2]
        center[j,i] = solp[1]
        print(i,j,p,solp)
        
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--fit-fwhm.fits', fwhm, hdr, clobber=True)
pyfits.writeto('bquint_datared/cal6600-20_013--phc-residual--fit-enter.fits', center, hdr, clobber=True)
get_ipython().magic('logstop')
