#!/usr/bin/python
#-*- codign: utf8 -*-

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
    
    parser.add_argument('-n', '--npoints', default=50, type=int, 
                        help="Number of points that will be used to fit" +
                        "the phase-map [50]")
    
    parser.add_argument('-o', '--output', type=str, default=None, 
                        help="Name of the output phase-map file.")
    
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Run program quietly.")
    
    args = parser.parse_args()
    v = not args.quiet
    if v:
        print("\n Phase-Map Fitting for BTFI")
        print(" by Bruno Quint & Fabricio Ferrari")
        print(" version 0.0a - Jan 2014")
        
    # Loading observed map ----------------------------------------------------
    if v: 
        print("\n Loading file: %s" % args.filename)        
    phase_map = pyfits.open(args.filename)[0]
    
    # Check if file was obtained with BTFI instrument
    header = phase_map.header
    if header['INSTRUME'].upper() not in ['BTFI'] and v:
        if v: print(" [Warning]: %s file was not obtained with BTFI instrument."
                     % args.filename)

    # Check whether data-cube was obtained from a FP scan or an iBTF scan
    if header['INSTRMOD'].upper() in ['IBTF'] and v:
        if v: print(" File obtained through an iBTF scan.")
    
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
        print("")
            
    # Fitting phase-map for a Fabry-Perot Map ---------------------------------
    elif header['INSTRMOD'].upper() in ['FP', 'FABRY-PEROT']:
        if v: 
            print(" File obtained through a Fabry-Perot scan.")
        
        width = header['NAXIS1']
        height = header['NAXIS2']
        ref_x = header['PHMREFX']
        ref_y = header['PHMREFY']
        unit = header['PHMUNIT']
        phmap = phase_map.data
        
        x = (numpy.linspace(0.10, 0.90, 50) * width).astype(int) 
        y = (numpy.linspace(0.10, 0.90, 50) * height).astype(int)
        X, Y = numpy.meshgrid(x, y)
        R = numpy.sqrt((X - ref_x) ** 2 + (Y - ref_y) ** 2)
        Z = phmap[Y, X]
        
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
        
        plt.figure()
        plt.title('Radial Plot of the Phase-Map')
        plt.plot(r, z, 'k.')
        plt.axhline((1.05 * z[0]), c='red')
        plt.xlabel('Radius [px]')
        plt.ylabel('Peak displacement [%s]' % unit)
        plt.grid()
        
        condition = numpy.where(z <= (1.05 * z[0]), True, False)
        not_condition = numpy.where(z <= (1.05 * z[0]), False, True)
        not_r = r[not_condition]
        not_z = z[not_condition]
        r = r[condition]
        z = z[condition]
                
        fit_func = lambda P, R: P[0] + P[1] * R + P[2] * R ** 2
        err_func = lambda P, R, Z: Z - fit_func(P, R)
        p, _ = optimize.leastsq(err_func, [0, 0, 0], args=(r, z))
        rr = numpy.linspace(r.min(), r.max(), 1000)
        zz = fit_func(p, rr)
        
        FSR = numpy.median(numpy.abs(not_z - fit_func(p, not_r)))
        print(" FSR = %f" % FSR)
                
        plt.figure()
        plt.title('Radial Plot of the Phase-Map')
        plt.plot(r, z, '.', color="#0000FF")
        plt.plot(not_r, not_z - FSR, '.', color="#9999FF")
        plt.plot(rr, zz, 'r-', lw=2)
        plt.xlabel('Radius [px]')
        plt.ylabel('Peak displacement [%s]' % unit)
        plt.grid()
    
        e = z - fit_func(p, r) 
        if v:
            print("  phi(x,y) = %.2e + %.2e x + %.2e x^2" % (p[0], p[1], p[2])) 
            print("  Error abs min: %f" % numpy.abs(e).min())
            print("  Error avg: %f" % e.mean())
            print("  Error std: %f" % e.std())
            print("  Error rms: %f" % numpy.sqrt(((e ** 2).mean())))
            print("  Sampling in Z: %s" % phase_map.header['fpzdelt'])
                   
        header.set('PHMFITSR', FSR, 'Free spectral range')
        header.set('', '', before='PHMFITSR')
        header.set('', '--- Phase-Map Fitting ---', before='PHMFITSR')
        
        r, z = numpy.append(r, not_r), numpy.append(z, not_z - FSR)
        p, _ = optimize.leastsq(err_func, p, args=(r, z))
        rr = numpy.linspace(r.min(), r.max(), 1000)
        zz = fit_func(p, rr)
        
        plt.figure()
        plt.title('Radial Plot of the Phase-Map')
        plt.plot(r, z, 'k.')
        plt.plot(rr, zz, 'r-', lw=2)
        plt.xlabel('Radius [px]')
        plt.ylabel('Peak displacement [%s]' % unit)
        plt.grid()

        x = numpy.arange(width) 
        y = numpy.arange(height)
        X, Y = numpy.meshgrid(x, y)
        R = numpy.sqrt((X - ref_x) ** 2 + (Y - ref_y) ** 2)
        Z = fit_func(p, R)
        Z = Z - Z[ref_y, ref_x]
        phmap = phmap - phmap[ref_y, ref_x]
        
        fname = header['PHMREFF']
        fname = os.path.splitext(fname)[0]
        pyfits.writeto(fname + '--fit_phmap.fits', Z, header, clobber=True)
        pyfits.writeto(fname + '--res_phmap.fits', Z - phmap, header, clobber=True)
          
        plt.show() 
        
    else:
        if v: print(" [Warning]: File was not obtained from FP or iBTF.")
        if v: print(" [Warning]: Don't know what to do. Leaving now.\n")
        from sys import exit
        exit()
    
#==============================================================================
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

if __name__ == '__main__':
    main()
