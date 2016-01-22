"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import astropy.io.fits as pyfits

FFMpegWriter = animation.writers['ffmpeg']

filename = "/home/bquint/Dropbox/new_cube30Dor_SOAR_350channels.fits"
data = pyfits.getdata(filename)

fig, ax = plt.subplots()
im = plt.imshow(data[0, :, :], origin='lower', interpolation='nearest',
                cmap=plt.get_cmap('cubehelix_r'))

ax.grid()
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

def animate(i):
    im.set_array(data[i, :, :])
    return im,

ani = animation.FuncAnimation(fig, animate, np.arange(1, 350),
                              interval=25, blit=True)
plt.show()