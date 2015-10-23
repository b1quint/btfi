#!/usr/bin/python
#-*- codign: utf8 -*-
"""
    2014.04.16 15:51 - Fixed keyword to access phase-map sampling.
    2014.08.25 18:30 - Fixed ref_x/ref_y
                     - Show phmap, ref_x/ref_y, npoints
"""
from __future__ import division, print_function

def main():

    import argparse
    import matplotlib.pyplot as plt
    import numpy
    import os
    import scipy.optimize as optimize
    from astropy.io import fits as pyfits

    # Parsing arguments -------------------------------------------------------
    parser = argparse.ArgumentParser(description="Fits an existing phase-map.")

    parser.add_argument('filename',
                        type=str,
                        help="Input phase-map name.")

    parser.add_argument('-i', '--interactions', default=5, type=int,
                        help="Number of interactions in the process [5]")

    parser.add_argument('-n', '--npoints', default=2500, type=int,
                        help="Number of points that will be used to fit" +
                        "the phase-map [50]")

    parser.add_argument('-o', '--output', type=str, default=None,
                        help="Name of the output phase-map file.")

    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run program quietly.")

    parser.add_argument('-s', '--show_plots', action='store_true',
                        help="Show plots (good for checking quality of the observed phase-map and the fitting.")

    args = parser.parse_args()
    v = not args.quiet
    if v:
        print("\n Phase-Map Fitting for BTFI")
        print(" by Bruno Quint & Fabricio Ferrari")
        print(" version 0.0a - Jan 2014")
        print("")

    check_dimensions(args.filename, dimensions=2)

    # Loading observed map ----------------------------------------------------
    if v:
        print(" Loading file: %s" % args.filename)
    phase_map = pyfits.open(args.filename)[0]

    # Check if file was obtained with BTFI instrument
    header = phase_map.header
    try:
        if header['INSTRUME'].upper() not in ['BTFI'] and v:
            if v: print(" [Warning]: %s file was not obtained with BTFI instrument."
                         % args.filename)
    except KeyError:
        warning("'INSTRUME' card was not found in the files' header.")

    mode = check_mode(args.filename)

    # Fitting Phase-Map for IBTF ----------------------------------------------
    if mode == 'ibtf':
        if v:
            print(" File obtained through an iBTF scan.")

        width = phase_map.header['naxis1']
        height = phase_map.header['naxis2']
        vmin = phase_map.data.mean() - 1.5 * phase_map.data.std()
        vmax = phase_map.data.mean() + 1.5 * phase_map.data.std()
        plt_config = {'origin': 'lower',
                      'cmap': get_colormap(),
                      'interpolation': 'nearest',
                      'vmin': vmin, 'vmax': vmax}
        if v:
            print(" Phase-map dimensions: [%d, %d]" % (width, height))
            print(" Done.")

        plt.subplot(131)

        plt.imshow(phase_map.data, **plt_config)
        plt.xticks([]), plt.yticks([])
        plt.xlabel('Observed Map')

        # Starting fitting process ------------------------------------------------
        npoints = numpy.sqrt(args.npoints).astype(int)
        if v:
            print("\n Starting phase-map fitting.")
            print(" %d x %d points will be used in the process." %
                   (npoints, npoints))

        x = (numpy.linspace(0.1, 0.9, npoints) * width).astype(int)
        y = (numpy.linspace(0.1, 0.9, npoints) * height).astype(int)
        x, y = numpy.meshgrid(x, y)
        x = numpy.ravel(x)
        y = numpy.ravel(y)
        z = numpy.ravel(phase_map.data[y,x])

        fit_func = lambda p, x, y: p[0] + p[1] * x + p[2] * y
        err_func = lambda p, x, y, z: z - fit_func(p, x, y)
        params = [z.mean(), 0, 0]

        # Fitting Plane -----------------------------------------------------------
        X = numpy.arange(phase_map.header['naxis1'])
        Y = numpy.arange(phase_map.header['naxis2'])
        X, Y = numpy.meshgrid(X, Y)

        if v: print("")
        for i in range(args.interactions):
            if v:
                print(" Fitting plane - Interaction %d" % (i + 1))

            if i == 0: e = z
            condition = numpy.where(numpy.abs(e - e.mean()) <= e.std())
            xx = x[condition]
            yy = y[condition]
            zz = z[condition]

            params, _ = optimize.leastsq(err_func, params, args=(xx, yy, zz))

            Z = fit_func(params, X, Y)
            error = Z - phase_map.data
            e = numpy.ravel(error[y,x])

        if v:
            p = params
            print("  phi(x,y) = %.2f + %.2fx + %.2fy" % (p[0], p[1], p[2]))
            print("  Error abs min: %f" % numpy.abs(e).min())
            print("  Error avg: %f" % e.mean())
            print("  Error std: %f" % e.std())
            print("  Error rms: %f" % numpy.sqrt(((e ** 2).mean())))

        plt.scatter(xx, yy, c=zz, cmap=get_colormap())
        plt.xlim(0, width), plt.ylim(0, height)

        plt.subplot(132)
        plt.imshow(error, **plt_config)
        plt.xticks([]), plt.yticks([])
        plt.xlabel("Residual")

        plt.subplot(133)
        plt.imshow(Z, **plt_config)
        plt.xticks([]), plt.yticks([])
        plt.xlabel("Fitted map")

        plt.show()

        ref_x = header['PHMREFX']
        ref_y = header['PHMREFY']
        fname = header['PHMREFF']

        fname = os.path.splitext(fname)[0]
        pyfits.writeto(fname + '--fit_phmap.fits', Z - Z[ref_y, ref_x], header, clobber=True)
        pyfits.writeto(fname + '--res_phmap.fits', Z - phase_map.data, header, clobber=True)

        print("")

    # Fitting phase-map for a Fabry-Perot Map ---------------------------------
    elif mode == 'fp':
        npoints = numpy.sqrt(args.npoints).astype(int)

        if v:
            print(" File obtained through a Fabry-Perot scan.")
            print(" Starting phase-map fitting.")
            print(" %d x %d points will be used in the process." %
                   (npoints, npoints))

        # Read data
        width = header['NAXIS1']
        height = header['NAXIS2']
        ref_x = header['PHMREFX']
        ref_y = header['PHMREFY']
        unit = header['PHMUNIT']
        sampling = header['PHMSAMP']
        FSR = header['PHMFSR']
        phmap = phase_map.data

        # From coordinates to pixels
        try:
            ref_x = header['CRPIX1'] + (ref_x - header['CRVAL1']) / header['CDELT1']
            ref_y = header['CRPIX2'] + (ref_y - header['CRVAL2']) / header['CDELT2']
        except KeyError:
            warning(" WCS not found. Using phisical coordinates.")

        phmap = phmap - phmap[ref_y, ref_x]

        x = (numpy.linspace(0.05, 0.95, npoints) * width).astype(int)
        y = (numpy.linspace(0.05, 0.95, npoints) * height).astype(int)
        X, Y = numpy.meshgrid(x, y)
        R = numpy.sqrt((X - ref_x) ** 2 + (Y - ref_y) ** 2)
        Z = phmap[Y, X]

        if args.show_plots:
            phmap_figure = plt.figure()
            phmap_axes = phmap_figure.add_subplot(111)
            phmap_imshow = phmap_axes.imshow(phmap, origin='lower',
                                             interpolation='nearest',
                                             cmap='coolwarm')
            phmap_axes.scatter(ref_x, ref_y, c='orange', s=400, marker="+", lw=3)
            phmap_axes.scatter(X, Y, c='g', s=1, marker=".", alpha=0.7)
            phmap_axes.set_xlabel("X [px]")
            phmap_axes.set_ylabel("Y [px]")
            phmap_axes.set_xlim(0, width)
            phmap_axes.set_ylim(0, height)
            phmap_axes.grid()
            phmap_figure.colorbar(phmap_imshow)


        x = numpy.ravel(X)
        y = numpy.ravel(Y)
        r = numpy.sqrt((x - ref_x) ** 2 + (y - ref_y) ** 2)
        z = numpy.ravel(Z)

        condition = numpy.where(z > z.min(), True, False) * \
                    numpy.where(z < z.max(), True, False)

        r = r[condition]
        z = z[condition]

        z = z[numpy.argsort(r)]
        r = numpy.sort(r)

        # Checking if parabola is up or down.
        if v:
            print("\n Checking if parabola is up or down.")

        dz = numpy.diff(z,1)
        dz_abs = numpy.abs(dz)
        dz_sign = numpy.sign(dz)
        sign = numpy.median(dz_sign[(dz_sign != 0) * (dz_abs <= sampling)])

        if v:
            print("  Parabola is %s" % ('up' if sign > 0 else 'down'))

        # Tell me the limits to fit the first parabola
        where = numpy.argmin(numpy.abs(r[dz_abs >= FSR / 2][0] - r))

        # Plot the gradient
        if args.show_plots:
            plt.figure(figsize=(16,7))
            plt.subplot(2,2,3)
            plt.plot(r[1:], dz, 'b-')
            plt.gca().yaxis.set_label_position("right")
            plt.axvline(r[where], color='black', lw=2, ls='--')
            plt.axhline(FSR / 2, color='red', ls='--', label="FSR")
            plt.axhline(- FSR / 2, color='red', ls='--')
            plt.xlabel('Radius [px]')
            plt.ylabel('Gradient \n [%s]' % unit)
            plt.legend(loc='best')
            plt.grid()

        # This is the first fit
        p = numpy.polyfit(r[:where], z[:where], 2)
        rr = numpy.linspace(r[0], r[where], 1000)
        zz = numpy.polyval(p, rr)

        # Plot the data before correction
        if args.show_plots:
            plt.subplot(2,2,1)
            plt.plot(r[:where], z[:where], 'b.', alpha=0.25, label='Not to be fixed')
            plt.plot(r[where:], z[where:], 'r.', alpha=0.25, label='Data to be fixed')
            plt.plot(rr, zz, 'k-', lw=2)
            plt.axvline(r[where], color='black', lw=2, ls='--')
            plt.gca().yaxis.set_label_position("right")
            plt.xlabel('Radius [px]')
            plt.ylabel('Peak displacement \n [%s]' % unit)
            plt.legend(loc='best')
            plt.grid()

        # Displace the FSR
        error = numpy.abs(z - numpy.polyval(p, r) + sign * FSR)

        # Plot error
        if args.show_plots:
            plt.subplot(2,2,4)
            plt.plot(r, error, 'k.', alpha=0.25)
            # plt.gca().yaxis.tick_right()
            plt.gca().yaxis.set_label_position("right")
            plt.xlabel('Radius [px]')
            plt.ylabel('Error \n [%s]' % unit)
            plt.ylim(ymin=-50, ymax=1.1*error.max())
            plt.grid()

        condition = (error > 2 * sampling)

        # Plot data after correction
        if args.show_plots:
            plt.subplot(2,2,2)
            plt.plot(r[condition], z[condition], 'b.', alpha=0.25,
                     label='Not fixed data')
            plt.plot(r[~condition], z[~condition] + sign * FSR, 'r.',
                     alpha=0.25, label='Fixed data')
            plt.gca().yaxis.set_label_position("right")
            plt.xlabel('Radius [px]')
            plt.ylabel('Peak displacement \n [%s]' % unit)
            plt.grid()

        # This is the second fit
        z = numpy.where(error >= 2 * sampling, z, z + sign * FSR)
        p = numpy.polyfit(r, z, 2)

        if args.show_plots:
            rr = numpy.linspace(r[0], r[-1], 1000)
            zz = numpy.polyval(p, rr)
            plt.plot(rr, zz, 'k-', lw=2, label='Fitted data.')
            plt.legend(loc='best')

        error = z - numpy.polyval(p, r)

        if v:
            print("  phi(x,y) = %.2e x^2 + %.2e x + %.2e " % (p[0], p[1], p[2]))
            print("  Error abs min: %f" % numpy.abs(error).min())
            print("  Error avg: %f" % error.mean())
            print("  Error std: %f" % error.std())
            print("  Error rms: %f" % numpy.sqrt(((error** 2).mean())))
            print("  Sampling in Z: %s" % phase_map.header['phmsamp'])
            print(" ")

        x = numpy.arange(width)
        y = numpy.arange(height)
        X, Y = numpy.meshgrid(x, y)
        R = numpy.sqrt((X - ref_x) ** 2 + (Y - ref_y) ** 2)
        Z = numpy.polyval(p, R)
        Z = Z - Z[ref_y, ref_x]

        fname = header['PHMREFF']
        fname = os.path.splitext(fname)[0]
        pyfits.writeto(fname + '--fit_phmap.fits', Z, header, clobber=True)
        pyfits.writeto(fname + '--res_phmap.fits', Z - phmap, header, clobber=True)

        if v:
            print(" All done.\n")

        if args.show_plots:
            plt.show()

    else:
        if v: print(" [Warning]: File was not obtained from FP or iBTF.")
        if v: print(" [Warning]: Don't know what to do. Leaving now.\n")
        from sys import exit
        exit()

def check_dimensions(filename, dimensions=3, keyword='NAXIS'):
    """
    Method written to check the dimensions of the input fits data.
    """
    from astropy.io import fits as pyfits
    import sys

    header = pyfits.getheader(filename)

    if not keyword in header:
        data = pyfits.getdata(filename)
        ndim = data.ndim
    else:
        ndim = header[keyword]

    if ndim is not dimensions:
        error(" INVALID DATA: wrong number of dimensions")
        error(" Leaving now!")
        print("")
        sys.exit()
    else:
        return

def check_mode(filename, keyword='INSTRMOD'):
    """
    Return if BTFI was obtained with a Fabry-Perot or with the iBTF.
    """
    from astropy.io import fits as pyfits

    header = pyfits.getheader(filename)

    if keyword not in header:
        warning("Instrument mode not found.")
        instrument_mode = ''
        while instrument_mode.lower() not in ['ibtf', 'fp']:
            instrument_mode = raw_input(" Enter 'ibtf' or 'fp': ")
    else:
        if header[keyword].upper() in ['IBTF']:
            instrument_mode = 'ibtf'

        if header[keyword].upper() in ['FP', 'FABRY-PEROT']:
            instrument_mode = 'fp'

    return instrument_mode

def get_colormap():

    from matplotlib import colors

    cdict = {'red': ((0.0, 1.0, 1.0),
                     (0.25, 1.0, 1.0),
                     (0.5, 0.0, 0.0),
                     (0.75, 0.0, 0.0),
                     (1.0, 0.9, 0.9)),
            'green': ((0.0, 0.9, 0.9),
                     (0.25, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (0.75, 0.0, 0.0),
                     (1.0, 0.9, 0.9)),
            'blue': ((0.0, 0.9, 0.9),
                     (0.25, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (0.75, 1.0, 1.0),
                     (1.0, 1.0, 1.0))}

    return colors.LinearSegmentedColormap('heaven_hell',cdict,256)

def error(my_string):
    s = bcolors.FAIL + '[ERROR] ' + bcolors.ENDC
    s = s + str(my_string)
    print(s)
    return

def warning(my_string):
    s = bcolors.WARNING + '[!] ' + bcolors.ENDC
    s = s + str(my_string)
    print(s)
    return

class bcolors:

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

if __name__ == '__main__':
    main()
