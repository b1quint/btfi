#!/bin/csh -f

# General parameters:
#      - You requested to use the following FP: Thickness = 200 microns; tunable gap = 2 microns; p=609 @ Halpha 
#      - You requested to do a OBSERVATION (and not a calibration)
#      - The name of the object                         : NGC2440pose2
#      - The wavelength (at rest) you gave is             = 6562.78 angstroms
#      - The radial velocity is                           = 0 km/s
#      - The wavelength (redshifted)                      = 6562.78 angstroms
# Interference order:
#      - The interference order @ 6562.78                 = 609.498 
#      - The interference order @ 6598.95                 = 606.157 
#      - The interference order @ 6562.78                 = 609.498 
# Free Spectral Range :
#      - The FSR @ 6598.95 in wavelength                  = 10.8866 Angstrom
#      - The FSR @ 6562.78 in wavelength                  = 10.7675 Angstrom
#      - The FSR @ 6562.78 in thickness                   = 0.328139 microns 
#      - The FSR @ 6562.78 in wavelength                  = 10.7675 Angstrom
#      - The FSR @ 6562.78 in km/s                        = 491.869 km/s
#      - The queensgate constant QGC                      = 18.7072 Angstrom
#      - The FSR in BCV @ 6562.78A                        = 350.816
#      - The FSR in BCV @ 6562.78A                        = 350.816
#      - The FSR in BCV @ 6598.95A                        = 352.75
# Finesse & Scanning:
#      - You gave a real finesse                         = 17.75
#      - Shannon sampling of the finesse                 = 2.1
#      - Considering F=17.75 and the sampling =2.1, the float nb of ch to scan for one FSR  = 37.275
#      - Considering F=17.75 and FSR=10.7675, the spectral sampling = 0.606622 Angstroms
#      - The spectral Resolution @ 6562.78 Angstroms        = 10818
#      - The average number of BCV for one FSR             = 9.41157
# Overscanning:
#      - You wanted to scan 1 FSR 
#      - The BCV gap that will be scanned @ 6562.78 Angstro = 350.816
#      - The total number of channels that will be scanned  = 38
#      - The initial BCV value   (nfiniz)                   = 768
#      - The final BCV value should be around (nfiniz_end)  = 419.772
# SAMI:
#      - You gave nsweeps  = 1
#      - You gave nsteps   = 1
#      - You gave nframe   = 1
#      - You gave exptim per channel             = 10 seconds
#      - Readout time per exposure               = 2 seconds 
#      - Total exposure time (whole observation) = 7.6 minutes
#      - Total exposure time (whole observation) = 0.126667 hours
#      - You gave binxy                          = 4 
#      - You gave the basename                   = fp_obs

set dat = `date +%Y-%m-%dT%H:%M:%S`
set scid = "SCAN_$dat"
echo "SCAN $scid"
set sweepkey = "FAPERSWP"
set stepkey = "FAPERSST"
set scankey = "FAPERSID"
set nsweeps = 1
set nsteps = 1
set nframe = 1
set nfiniz = 768
set exptim = 10.0
set binxy = 4
set basename = "fp_obs"
set cmd = `sami dhe set image.dir /home2/images/20151204/015`
set cmd = `sami dhe dbs set $scankey $scid`
set cmd = `sami dhe dbs set $stepkey custom`
echo "setting number of images, exposure time and basename"
sami dhe set binning $binxy $binxy
sami dhe set obs.nimages $nframe
sami dhe set obs.exptime $exptim
sami dhe set image.basename $basename
echo
echo "image $basename, exptime $exptim"
echo "binning $binxy"
echo
echo "moving FP to channel 1: BCV=768"
sami FP moveabs 768
set sweepid = C001
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 1)"
sami dhe expose
echo
echo "moving FP to channel 2: BCV=759"
sami FP moveabs 759
set sweepid = C002
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 2)"
sami dhe expose
echo
echo "moving FP to channel 3: BCV=749"
sami FP moveabs 749
set sweepid = C003
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 3)"
sami dhe expose
echo
echo "moving FP to channel 4: BCV=740"
sami FP moveabs 740
set sweepid = C004
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 4)"
sami dhe expose
echo
echo "moving FP to channel 5: BCV=730"
sami FP moveabs 730
set sweepid = C005
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 5)"
sami dhe expose
echo
echo "moving FP to channel 6: BCV=721"
sami FP moveabs 721
set sweepid = C006
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 6)"
sami dhe expose
echo
echo "moving FP to channel 7: BCV=712"
sami FP moveabs 712
set sweepid = C007
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 7)"
sami dhe expose
echo
echo "moving FP to channel 8: BCV=702"
sami FP moveabs 702
set sweepid = C008
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 8)"
sami dhe expose
echo
echo "moving FP to channel 9: BCV=693"
sami FP moveabs 693
set sweepid = C009
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 9)"
sami dhe expose
echo
echo "moving FP to channel 10: BCV=683"
sami FP moveabs 683
set sweepid = C010
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 10)"
sami dhe expose
echo
echo "moving FP to channel 11: BCV=674"
sami FP moveabs 674
set sweepid = C011
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 11)"
sami dhe expose
echo
echo "moving FP to channel 12: BCV=664"
sami FP moveabs 664
set sweepid = C012
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 12)"
sami dhe expose
echo
echo "moving FP to channel 13: BCV=655"
sami FP moveabs 655
set sweepid = C013
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 13)"
sami dhe expose
echo
echo "moving FP to channel 14: BCV=646"
sami FP moveabs 646
set sweepid = C014
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 14)"
sami dhe expose
echo
echo "moving FP to channel 15: BCV=636"
sami FP moveabs 636
set sweepid = C015
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 15)"
sami dhe expose
echo
echo "moving FP to channel 16: BCV=627"
sami FP moveabs 627
set sweepid = C016
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 16)"
sami dhe expose
echo
echo "moving FP to channel 17: BCV=617"
sami FP moveabs 617
set sweepid = C017
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 17)"
sami dhe expose
echo
echo "moving FP to channel 18: BCV=608"
sami FP moveabs 608
set sweepid = C018
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 18)"
sami dhe expose
echo
echo "moving FP to channel 19: BCV=599"
sami FP moveabs 599
set sweepid = C019
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 19)"
sami dhe expose
echo
echo "moving FP to channel 20: BCV=589"
sami FP moveabs 589
set sweepid = C020
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 20)"
sami dhe expose
echo
echo "moving FP to channel 21: BCV=580"
sami FP moveabs 580
set sweepid = C021
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 21)"
sami dhe expose
echo
echo "moving FP to channel 22: BCV=570"
sami FP moveabs 570
set sweepid = C022
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 22)"
sami dhe expose
echo
echo "moving FP to channel 23: BCV=561"
sami FP moveabs 561
set sweepid = C023
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 23)"
sami dhe expose
echo
echo "moving FP to channel 24: BCV=552"
sami FP moveabs 552
set sweepid = C024
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 24)"
sami dhe expose
echo
echo "moving FP to channel 25: BCV=542"
sami FP moveabs 542
set sweepid = C025
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 25)"
sami dhe expose
echo
echo "moving FP to channel 26: BCV=533"
sami FP moveabs 533
set sweepid = C026
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 26)"
sami dhe expose
echo
echo "moving FP to channel 27: BCV=523"
sami FP moveabs 523
set sweepid = C027
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 27)"
sami dhe expose
echo
echo "moving FP to channel 28: BCV=514"
sami FP moveabs 514
set sweepid = C028
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 28)"
sami dhe expose
echo
echo "moving FP to channel 29: BCV=504"
sami FP moveabs 504
set sweepid = C029
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 29)"
sami dhe expose
echo
echo "moving FP to channel 30: BCV=495"
sami FP moveabs 495
set sweepid = C030
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 30)"
sami dhe expose
echo
echo "moving FP to channel 31: BCV=486"
sami FP moveabs 486
set sweepid = C031
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 31)"
sami dhe expose
echo
echo "moving FP to channel 32: BCV=476"
sami FP moveabs 476
set sweepid = C032
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 32)"
sami dhe expose
echo
echo "moving FP to channel 33: BCV=467"
sami FP moveabs 467
set sweepid = C033
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 33)"
sami dhe expose
echo
echo "moving FP to channel 34: BCV=457"
sami FP moveabs 457
set sweepid = C034
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 34)"
sami dhe expose
echo
echo "moving FP to channel 35: BCV=448"
sami FP moveabs 448
set sweepid = C035
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 35)"
sami dhe expose
echo
echo "moving FP to channel 36: BCV=439"
sami FP moveabs 439
set sweepid = C036
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 36)"
sami dhe expose
echo
echo "moving FP to channel 37: BCV=429"
sami FP moveabs 429
set sweepid = C037
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 37)"
sami dhe expose
echo
echo "moving FP to channel 38: BCV=420"
sami FP moveabs 420
set sweepid = C038
set cmd = `sami dhe dbs set $sweepkey $sweepid`
sami dhe set image.basename $basename"_"$sweepid
echo "SWEEP $sweepid"
echo "taking data...(sweep $sweepid step 38)"
sami dhe expose
# Channel: +Step ==> BCV
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]
# [0, -9, -10, -9, -10, -9, -9, -10, -9, -10, -9, -10, -9, -9, -10, -9, -10, -9, -9, -10, -9, -10, -9, -9, -10, -9, -10, -9, -10, -9, -9, -10, -9, -10, -9, -9, -10, -9]
# [768, 759, 749, 740, 730, 721, 712, 702, 693, 683, 674, 664, 655, 646, 636, 627, 617, 608, 599, 589, 580, 570, 561, 552, 542, 533, 523, 514, 504, 495, 486, 476, 467, 457, 448, 439, 429, 420]
