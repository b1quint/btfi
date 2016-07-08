#!/bin/env python
import time
import datetime
import sys
import os
"""
    
    THIS PROGRAM COMPUTE A SCANNING SEQUENCE FOR PF/SAM/SOAR
    Philippe Amram
    last update: 2016-04-02 By B. Quint

NOTATIONS:
    epais = distance between the two plates
    gap = the maximum tuning gap
    QGC = Queensgate Constant
    BCV = Binary Control Value
    
    INTERACTIVE
    In interactive mode, interactive = True, in non-interactive mode, interactive = False
"""

print("\n%s" % ("-"*100))
print("\n WELCOME ! ")
print(" %s" % time.strftime('%a, %d %b %Y %H:%M:%S GMT'))
print("\n     This program prepares your script to phasemap_fit on FP/SAMI")

interactive = True
#interactive = False

def main():

    """
   
CONSTANTS
   
    """
    celerite    = 299792.458 # Speed of light [km/s]
    lamb_halpha = 6562.78 # [A]
    lamb_Ne     = 6598.9529 # [A]
    bcv_max     = 4095    # 4096 value starting from 0 up to 4095

    """
    
INITIALISATION OF THE TWO SCRIPT FILES
    
    """

    """ 1) INITIALISATION OF THE SCANNING SCRIPT """

    tt = time.strftime('%Y-%m-%dT%Hh%Mm%Ss')
    #time.struct_time(tm_year=2015, tm_mon=1, tm_mday=28, tm_hour=19, tm_min=1, tm_sec=11, tm_wday=2, tm_yday=28, tm_isdst=0)

    """ 2) INITIALISATION OF THE RUNNING DIRECTORY """

    #dirtime = time.strftime('%Y%m%d')
    dirtime = "20160402"
    if interactive:
        print("\n Data of the observing phasemap_fit. ")
        print(" The date of the phasemap_fit (e.g. 20150130) will be automatically added to the address of the directory you will give now")
        sdir = raw_input(" Please, input the directory name (e.g.:/home2/images/): ")
        if sdir == "":
            sdir = '/home2/images/'
        sdir = os.path.abspath(sdir)

        print("     Your input is: %s" % sdir)
        print("     The data of the day will go into the directory : ", os.path.join(sdir, dirtime))
        running = raw_input("\n Give the running directory name where you will put and phasemap_fit the script (e.g. 001): ")
        running = str(running)
        sdir = os.path.join(sdir, dirtime, running)
        print("     The name of the directory where the script will be ran is : ",sdir)
    else:
        print(dirtime)
        running = input("\n Give the running directory name where you will put and phasemap_fit the script (e.g. 001): ")
        #running="012"
        sdir = os.path.join("/home2/images/", dirtime, running)

    """ 3) SCRIPT TO RUN TO COPY THE SCANNING SCRIPT FROM MY COMPUTER TO BTFI COMPUTER """

    tt = running
    ttsh = tt+'.sh'
    Fichier = open(ttsh,'w')
    Fichier.write("#!/bin/csh -f\n\n")
    Fichier0 = open('scpbtfidr.sh','w')
    Fichier0.write("#!/bin/csh -f\n\n")
    Fichier0.write("sshpass -p \"btfi88\" scp %s btfidr@139.229.18.227:/data/%s/%s/\n" % (ttsh,dirtime,running))
    # btfidr@btfidr:/data/20150317
    Fichier0.close()

    """
    
FABRY-PEROT TO USE
    
    """

    #dico0 = {'FP':["Thickness =  44 microns; tunable gap = 2 microns; p=134 @ Halpha","Thickness = 200 microns; tunable gap = 2 microns; p=609 @ Halpha"]}
    dico0 = {}
    dico0[1,1] = "Thickness =  44 microns; tunable gap = 2 microns; p=134 @ Halpha"
    dico0[1,2] = "Thickness = 200 microns; tunable gap = 2 microns; p=609 @ Halpha"

    if interactive:
        print("\n Please, input the name of the Fabry-Perot you wanted to use: ")
        print("      For TF (%s) put (1) " % dico0[1,1])
        print("      For PF (%s) put (2) " % dico0[1,2])
        pftf = int(raw_input(" Your choise : "))
    else:
        pftf = 2

    if pftf > 2:
        print("      Sorry, you input a value not allowed, choose between (1) and (2), please restart ")
        sys.exit(0)
    else:
        if pftf == 1:
            print("     Your input is : %s" % dico0[1,1])
        if pftf == 2:
            print("     Your input is : %s" % dico0[1,2])

#epais = float(input(" Please, input the thickness (not the gap) of the interferometer in $mu m$ (it should be a float, eg epais=350): "))
#epais = p_order * lamb*1e-4 /2

    if pftf == 1:
        epais = 44.
    if pftf == 2:
        epais = 200.

    """
    
CALIBRATION OR OBSERVATION
    
    """

    if interactive:
        calibration = ''
        while calibration not in ['0', '1']:
            calibration = raw_input("\n Please, input if you make a calibration (0) or an observation (1): ")
        calibration = int(calibration)
        #calibration = 0
        if calibration == 0:
            lamb = lamb_Ne # Ne I
            print("     The wavelength of the Neon line will be: %.3f" % lamb)
    else:
        calibration = ''
        while calibration not in ['0', '1']:
            calibration = raw_input("\n Please, input if you make a calibration (0) or an observation (1): ")
        calibration = int(calibration)
                              
    if interactive:
        if calibration == 0:
            print("     You requested a calibration.")
            lamb = lamb_Ne
        if calibration == 1:
            print("     You resqueted an observation.")
            #lamb = lamb_halpha
            print("\n You have to give the wavelength at rest of the line you will observe.")
            lamb = float(raw_input(" Please, input this wavelength (in Angstrom, it should be a float, e.g. Hapha = 6562.78) : "))
            #lamb = lamb_halpha

    if not interactive:
        if calibration == 0:
            lamb = lamb_Ne # Ne I
            #lamb = 5460.742 # Green Hg
            #lamb = 4358.343 # Indigo Hg
            #lamb = 6506.5281 # Ne I
            #lamb = 6532.8822 # Ne I
        if calibration == 1:
            #lamb = 6590.4
            #lamb = 6571.64
            lamb = lamb_halpha
            #lamb = 7000.

    if lamb < 0:
        print("      Sorry, you input a value not allowed because %.3f should be greater than 0, please restart " % lamb)
        sys.exit(0)

    lamb_rest=lamb

    if calibration == 0:
        vitesse = 0 # km/s

    if calibration == 1:
        if interactive:
            object_name = raw_input("\n Please, input the name of your object (e.g. NGC 7331): ")
            print("     Your input is: %s " % object_name)
            vitesse = float(raw_input("\n Please, input the radial velocity of the galaxy (in km/s): "))
            print("     Your input is: %.2f km/s" % vitesse)
        else:
            object_name = raw_input("\n Please, input the name of your object (e.g. NGC 7331): ")
            print("     Your input is: %s " % object_name)
            vitesse = float(raw_input("\n Please, input the radial velocity of the galaxy (in km/s): "))
            print("     Your input is: %.2f km/s" % vitesse)
            #object_name = "Cartwheel"
            #vitesse = 9050 # km/s

    """
    
INITIAL PARAMETER COMPUTATION
    
    """

    def ISL(ll,pp):
        #fsr_lamb = lamb/p_order
        isl_ll  = ll/pp*(1+1/(pp*pp))
        return isl_ll
    
    def P_ORDER(ee,ll):
        porder = 2. * ee * 1E+4 / ll
        return porder
        
    lamb = (vitesse / celerite + 1) * lamb_rest

    p_order         = P_ORDER(epais,lamb)
    p_order_halpha  = P_ORDER(epais,lamb_halpha)
    p_order_Ne      = P_ORDER(epais,lamb_Ne)

    p_order0        = int(p_order)
    e_fsr           = epais /p_order

    fsr_lamb        = ISL(lamb,p_order)
    fsr_lamb_Ne     = ISL(lamb_Ne,p_order_Ne)
    fsr_lamb_Ha     = ISL(lamb_halpha,p_order_halpha)

    fsr_kms         = celerite * fsr_lamb / lamb

    Fichier.write("# General parameters:\n")
    Fichier.write("#      - You requested to use the following FP: %s} \n" % dico0[1,pftf])
    if calibration == 0 :
        Fichier.write("#      - You requested to do a CALIBRATION (and not an observation on the sky)\n")
    if calibration == 1 :
        Fichier.write("#      - You requested to do a OBSERVATION (and not a calibration)\n")
        Fichier.write("#      - The name of the object                         : %s\n" % object_name)
    Fichier.write("#      - The wavelength (at rest) you gave is             = %.3f angstroms\n" % lamb_rest)
    if calibration == 1 :
        Fichier.write("#      - The radial velocity is                           = %.3f km/s\n" % vitesse)
    if calibration == 1 :
        Fichier.write("#      - The wavelength (redshifted)                      = %.3f angstroms\n" %lamb)
    Fichier.write("# Interference order:\n")
    Fichier.write("#      - The interference order @ %.3f                 = %.3f \n" % (lamb_halpha,p_order_halpha))
    Fichier.write("#      - The interference order @ %.3f                 = %.3f \n" % (lamb_Ne,p_order_Ne))
    Fichier.write("#      - The interference order @ %.3f                 = %.3f \n" % (lamb,p_order))
    Fichier.write("# Free Spectral Range :\n")
    Fichier.write("#      - The FSR @ %.3f in wavelength                  = %.3f Angstrom\n" % (lamb_Ne,fsr_lamb_Ne))
    Fichier.write("#      - The FSR @ %.3f in wavelength                  = %.3f Angstrom\n" % (lamb_halpha,fsr_lamb_Ha))
    Fichier.write("#      - The FSR @ %.3f in thickness                   = %.3f microns \n" % (lamb,e_fsr))
    Fichier.write("#      - The FSR @ %.3f in wavelength                  = %.3f Angstrom\n" % (lamb,fsr_lamb))
    Fichier.write("#      - The FSR @ %.3f in km/s                        = %.3f km/s\n" % (lamb,fsr_kms))

    """
    
QUEENSGATE CONSTANT
    
    """

    if interactive:
        print("\n (1) If you know it, you can use the Queensgate Constant already measured with the SAME CS100 AND the the SAME FP.")
        print(" (2) If you do not know, you must put the total plate gap in BCV corresponding to one FSR at the wavelength")
        print("     '1' means you DO     want to give a the Queensgate Constant.")
        print("     '2' means you DO NOT want to give a the Queensgate Constant but a number of BCV corresponding to one FSR")
        QGC_or_not = int(input(" your input ('1' or '2'): "))
    else:
        QGC_or_not = 2

    if QGC_or_not > 2 or QGC_or_not < 1:
        print("      Sorry, you input %s which is a value not allowed please choose '1' or '2' " % (QG_or_not))
        sys.exit(0)

    if QGC_or_not == 1:
        if interactive:
            QGC = float(input("\n Please, input the Queensgate Constant (in Angstrom, could be a float, e.g. 9.30): "))
            #QGC = 9.15
            print("     Your input is: %s Angstroms." % (QGC))
        else:
            QGC = 9.15  # undersampling
            QGC = 9.40  # oversampling
            QGC = 9.30  # close to be perfect

    if QGC_or_not == 2:
        if interactive:
            print("\n You first must choose the wavelength at which the gap in BCV will be given.")
            print(" NOTE: this value is not necessary the scanning wavelength.")
            # lamb_QGC = float(input(" Give this wavelength (could be a float e.g. Ne 6598.9529): "))
            lamb_QGC = 6598.9529
            print("\n Please, input the total plate gap in BCV corresponding to one FSR at the wavelength %s of reference" % (lamb_QGC))
            print(" The BCV range between 0 and %s" % (bcv_max))
            fsr_bcv_lamb_QGC = float(input(" your input (could be a float, e.g. 705): "))
            #fsr_bcv_lamb_QGC = 352.75
            print("     Your input is: %s" % (fsr_bcv_lamb_QGC))
        else:
            fsr_bcv_lamb_QGC = 352.75
            lamb_QGC = lamb_Ne
        QGC = lamb_QGC / fsr_bcv_lamb_QGC
        print("     A queensgate has been computed : %s" % (QGC))

    fsr_bcv_lamb    = lamb / QGC
    fsr_bcv_lamb_Ha = lamb_halpha / QGC
    fsr_bcv_lamb_Ne = lamb_Ne / QGC

    """
    
NUMBER OF CHANNELS TO SCAN
    
    """

    if interactive:
        print("\n Taking into account the Finesse and the sampling, the number of channel to scan could be computed automatically.")
        print(" Alternatively you can define yourself the number of channels to scan.")
        print("     (1) You DO WISH to compute automatically the number of channels to scan")
        print("     (2) You DO NOT WISH to give manually the number of channels to scan")
        nchan_manuel = int(input("          Please give you choose (1 or 2): "))
    else:
        nchan_manuel = 1

    if nchan_manuel == 1:
        if interactive:
            finesse = float(input("\n Please, input the Finesse (finesse must be a float): "))
        else:
            finesse = 17.75
        if finesse <= 1:
            print("      Sorry, you input a value not allowed because %.3f should be greater than 1, please restart " % (finesse))
            sys.exit(0)
        if interactive:
            sampling = float(input("\n Please, input the sampling, Shannon indicates that the sampling could be 2 (could be a float): "))
        else:
            sampling = 2.1
        if (sampling) <= 1:
            print("      Sorry, you input a value not allowed because %.3f should be greater or equal to one, please restart " % (sampling))
            sys.exit(0)
        """ Integer value + 1 to avoid undersampling """
        nchan = sampling*finesse

    if nchan_manuel == 2:
        if interactive:
            nchan = int(input("\n Please input the number of channel to scan one FSR (must be an integer): "))
        else:
            nchan = 38

#    bcv_step = fsr_bcv_lamb / nchan / 2
    bcv_step = fsr_bcv_lamb / nchan
    if (bcv_step) < 2:
        print("\n     Sorry, your scanning step in BCV =%.3f is too small, it should not be lower than 2." % (bcv_step))
        if nchan_manuel == 1:
            print("     This could be due to the finesse (=%.3f) or/and the sampling (=%.3f) too high." % (finesse,sampling))
        if nchan_manuel == 2:
            print("     This could be due to the number of channels (=%.3f) too high." % (nchan))
        print("     Please RESTART from the beginning.")
        sys.exit(0)

    Fichier.write("#      - The queensgate constant QGC                      = %.3f Angstrom\n" % (QGC))
    Fichier.write("#      - The FSR in BCV @ %.3fA                        = %.3f\n" % (lamb,fsr_bcv_lamb))
    Fichier.write("#      - The FSR in BCV @ %.3fA                        = %.3f\n" % (lamb_halpha,fsr_bcv_lamb_Ha))
    Fichier.write("#      - The FSR in BCV @ %.3fA                        = %.3f\n" % (lamb_Ne,fsr_bcv_lamb_Ne))
    Fichier.write("# Finesse & Scanning:\n")
    if nchan_manuel == 1:
        Fichier.write("#      - You gave a real finesse                         = %.3f\n" % (finesse))
        Fichier.write("#      - Shannon sampling of the finesse                 = %.3f\n" % (sampling))
        Fichier.write("#      - Considering F=%.3f and the sampling =%.3f, the float nb of ch to scan for one FSR  = %.3f\n" % (finesse,sampling,nchan))
        Fichier.write("#      - Considering F=%.3f and FSR=%.3f, the spectral sampling = %.3f Angstroms\n" % (finesse,fsr_lamb,fsr_lamb/finesse))
        Fichier.write("#      - The spectral Resolution @ %.3f Angstroms        = %.3f\n" % (lamb,int(lamb*finesse/fsr_lamb)))
    else:
        Fichier.write("#      - The number of channels to scan for one FSR      = %.3f\n" % (nchan))
    Fichier.write("#      - The average number of BCV for one channel             = %.3f\n" % (bcv_step))

    """
    
SCAN MORE THAN ONE FSR ?
    
    """

    if interactive:
        print("\n You can scan more than one FSR.")
        print(" NOTE: The number of channel to scan for more than one FSR will be larger and computed automatically.")
        overlap = float(input(" Please, input the number of FSR you want to scan (could be a float, \"1\" means you will scan one FSR): "))
    else:
        overlap = 1.

    if overlap < 0:
        print("      Sorry, you input a value not allowed because %.3f should be greater than 0, please restart " % (overlap))
        sys.exit(0)

    if (fsr_bcv_lamb*overlap) > bcv_max:
        print("     \nSorry, you input a value not allowed because %.3f X %.3f = %.3f is greater than %.3f." % (int(fsr_bcv_lamb,overlap),int(fsr_bcv_lamb*overlap),bcv_max))
        print("     Please RESTART from the beginning.")
        sys.exit(0)
    else:
        fsr_bcv_lamb = fsr_bcv_lamb * overlap
        nchan = int(nchan * overlap)+1

    Fichier.write("# Overscanning:\n")
    Fichier.write("#      - You wanted to scan %.3f FSR \n" % (overlap))
    Fichier.write("#      - The BCV gap that will be scanned @ %.3f Angstro = %.3f\n" % (lamb,fsr_bcv_lamb))
    Fichier.write("#      - The total number of channels that will be scanned  = %.3f\n" % (nchan))

    """ TO SCAN IN DECREASING  THE RADIUS OF THE RINGS """
    #nfiniz0 = int(input(" Please, input the zero Z value (nfiniz0 must be an integer): "))
    #input(" Please, input the initial Z value (nfiniz must be an integer): "))
    #nfiniz0 =  0
    #       DIVIDED BY 4 BECAUSE OF THE UNAVAILABLE BCV RANGE
    #nfiniz0 =  int(bcv_max/4)
    #nfiniz = nfiniz0 - int(fsr_bcv_lamb/4.)

    """ TO SCAN IN INCREASING  THE RADIUS OF THE RINGS """
    #nfiniz0 =  int(bcv_max/4)
    #nfiniz = nfiniz0 + int(fsr_bcv_lamb/4.)
    nfiniz = 768
    nfiniz_end = nfiniz - (nchan - 1) * bcv_step

    """ Checking using the basic formula """
    base = lamb / QGC
    step = base / nchan
    # print("lamb= ",lamb," QGC =",QGC," nchan =",nchan," base (BCV)= ",base," step (BCV)= ",step)

    #Fichier.write("#      - The zero BCV value     (nfiniz0)                   = %.3f\n" % (nfiniz0))
    Fichier.write("#      - The initial BCV value   (nfiniz)                   = %.3f\n" % (nfiniz))
    Fichier.write("#      - The final BCV value should be around (nfiniz_end)  = %.3f\n" % (nfiniz_end))
    
    uneminute = 60. # second
    if (calibration == 0):
        basename = "fp_cal"
    if (calibration == 1):
        basename = "fp_obs"
    if interactive:
        #nsweeps    = int(input(" Please, input how many \"sweeps\" will be done on this scan (nsweeps must be an integer): "))
        nsweeps     = 1
        #nsteps     = int(input(" Please, input how many Z steps each sweep will have (nsteps must be an integer): "))
        nsteps      = 1
        #nframe     = int(input(" Please, input how many images we will take in each step (each Z value, nframe must be an integer): "))
        nframe      = 1
        basename    = raw_input("\n Please, set the basename of your fits image (basename must be a string, default: fp_sami): ")
        if basename == "":
            basename = "fp_sami"

        print("     Your basename is : %s" % basename)
        binxy       = input("\n Please, set the binning of the CCD image\n" + \
                            " (binxy must be an integer, e.g. 4 for a 4x4 binning)\n " + \
                            " Default 4\n > ")
        binxy = 4
        exptim      = float(input("\n Please, set the image exposure time per channel" + \
                                  " in seconds (exptim could be a float):\n > "))
        #exptim = 5
    else:
        nsweeps = 1
        binxy = 4
        nsteps = 1
        nframe = 1
        #basename = "fp_sami"
        exptim_min = 5
        exptim = exptim_min * uneminute
        exptim      = float(input("\n Please, set the image exposure time per channel" + \
                                  " in seconds (exptim could be a float):\n > "))

    readout_time = 2.        # 2 seconds = readout time @ binxy = 4 x 4 ???
    exptim_total = (nchan * (exptim + readout_time)) / uneminute

    if (exptim) < 0:
        print("      Sorry, you input a value not allowed because" + \
              " %.3f should be greater than 0, please restart " % (exptim))
        sys.exit(0)

    Fichier.write("# SAMI:\n")
    Fichier.write("#      - You gave nsweeps  = %s\n" % (nsweeps))
    Fichier.write("#      - You gave nsteps   = %s\n" % (nsteps))
    Fichier.write("#      - You gave nframe   = %s\n" % (nframe))
    Fichier.write("#      - You gave exptim per channel             = %.3f seconds\n" % (exptim))
    Fichier.write("#      - Readout time per exposure               = %.3f seconds \n" % (readout_time))
    Fichier.write("#      - Total exposure time (whole observation) = %.3f minutes\n" % (exptim_total))
    Fichier.write("#      - Total exposure time (whole observation) = %.3f hours\n" % (exptim_total/uneminute))
    Fichier.write("#      - You gave binxy                          = %s \n" % (binxy))
    Fichier.write("#      - You gave the basename                   = %s\n\n" % (basename))
    Fichier.write("set dat = `date +%Y-%m-%dT%H:%M:%S`\n")
    Fichier.write("set scid = \"SCAN_$dat\"\n")
    Fichier.write("echo \"SCAN $scid\"\n")
    Fichier.write("set sweepkey = \"FAPERSWP\"\n")
    Fichier.write("set stepkey = \"FAPERSST\"\n")
    Fichier.write("set scankey = \"FAPERSID\"\n")
    Fichier.write("set nsweeps = %s\n" % (nsweeps))
    Fichier.write("set nsteps = %s\n" % (nsteps))
    Fichier.write("set nframe = %s\n" % (nframe))
    Fichier.write("set nfiniz = %s\n" % (nfiniz))
    Fichier.write("set exptim = %s\n" % (exptim))
    Fichier.write("set binxy = %s\n" % (binxy))
    #Fichier.write("set basename = \"fp_sami\"\n")
    Fichier.write("set basename = \"%s\"\n" % (basename))
    Fichier.write("set cmd = `sami dhe set image.dir %s`\n" % (sdir))
    Fichier.write("set cmd = `sami dhe dbs set $scankey $scid`\n")
    Fichier.write("set cmd = `sami dhe dbs set $stepkey custom`\n")
    Fichier.write("echo \"setting number of images, exposure time and basename\"\n")
    Fichier.write("sami dhe set binning $binxy $binxy\n")
    Fichier.write("sami dhe set obs.nimages $nframe\n")
    Fichier.write("sami dhe set obs.exptime $exptim\n")
    Fichier.write("sami dhe set image.basename $basename\n")
    Fichier.write("echo\n")
    Fichier.write("echo \"image $basename, exptime $exptim\"\n")
    Fichier.write("echo \"binning $binxy\"\n")

    dico = {'channel':[], 'step':[], 'BCV':[]}
    iBCV = 0
    delta_iBCV = 0
    for cnt in range(1,nchan+1):
        iBCV0=iBCV
        #BCV = nfiniz + (cnt-1) * bcv_step
        BCV = nfiniz - (cnt-1) * bcv_step
        if BCV >= 0:
            if (int(BCV + 0.5) > int(BCV)):
                iBCV = int(BCV)+1
            else:
                iBCV = int(BCV)
        else:
            if (int(BCV - 0.5) < int(BCV)):
                iBCV = int(BCV)-1
            else:
                iBCV = int(BCV)
#       print("BCV=",BCV,"  ibcv=",iBCV)
        Fichier.write("echo\n")
        Fichier.write("echo \"moving FP to channel %s: BCV=%s\"\n" % (cnt,iBCV))
        Fichier.write("sami FP moveabs %s\n" % (iBCV))
#       Fichier.write("# At channel %s, BCV = %s, iBCV = %s" % (cnt,BCV,iBCV))
#       Fichier.write("sami fp moveabs $fpiniz\n")
#       Fichier.write("set sweepid = \"S\"$scnt\n")
        Fichier.write("set sweepid = C%03d\n"%cnt)
        Fichier.write("set cmd = `sami dhe dbs set $sweepkey $sweepid`\n")
        Fichier.write("sami dhe set image.basename $basename\"_\"$sweepid\n")
        Fichier.write("echo \"SWEEP $sweepid\"\n")
        Fichier.write("echo \"taking data...(sweep $sweepid step %s)\"\n" % (cnt))
        Fichier.write("sami dhe expose\n")
        if cnt == 1 :
            delta_iBCV = 0
        else:
            delta_iBCV = iBCV-iBCV0
        dico['channel'].append(cnt)
        dico['step'].append(delta_iBCV)
        dico['BCV'].append(iBCV)
    Fichier.write("# Channel: +Step ==> BCV\n")
    Fichier.write("# %s\n" % (dico['channel']))
    Fichier.write("# %s\n" % (dico['step']))
    Fichier.write("# %s\n" % (dico['BCV']))

    Fichier.close()

    print("\n     The name of the script you have to phasemap_fit on SAMI computer is : %s " % ttsh)
    print("     Copy the following script to SAMI computer in the following directory : %s " % sdir)
    print("     NOTE: You have to pass by BTFIDR computer to have access to SAMI computer")
    print("     To copy the script from your computer to BTFI computer,")
    print("     phasemap_fit the script \"scpbtfidr.sh\" which have been created now.")
    print("\n END OF THE PROGRAM")
    print(" %s" % ("-"*79))

if __name__ == '__main__':
    main()

