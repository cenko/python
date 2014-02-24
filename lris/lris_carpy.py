
# Global python packages
import pyraf
from pyraf import iraf
import os, shutil, math, subprocess
import pyfits

# Local python packages
from iqutils import *

# IRAF modules
iraf.images()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.specred()
iraf.rvsao()

# IRAF variables
yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets
imcombine=iraf.imcombine

# PyRAF setup
pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

# Global variables
LGRATINGS = {"400/8500": {"dispersion": 1.16, "trimreg": "*,*",
                          "slits": [2, 75, 510, 581, 1110]}}
LGRISMS = {"600/4000": {"dispersion": 0.63, "trimreg": "1421:2660,*",
                        "slits": [2, 41, 600, 641, 1200],
                        "wavelengths": [3010., 5600.]},
           "400/3400": {"dispersion": 1.09, "trimreg": "1421:2660,1451:4096",
                        "slits": [2, 41, 600, 641, 1200],
                        "wavelengths": [2850., 5740.]}}
LRARCLIST = "/Users/cenko/iraf/linelists/lris_carpy_listR.dat"
LBARCLIST = "/Users/cenko/iraf/linelists/lris_carpy_listB.dat"
LSKYLIST = "/Users/cenko/iraf/linelists/skylines.lris.txt"
RCCDRNOISE = 4.8 # CCD read noise (e-)
RCCDGAIN = 1.2 # CCD e- per DN
BCCDRNOISE = 4.2
BCCDGAIN = 1.6
BFLAT = "Flat-B.fits" # Median combined flat field for blue chip
BSHARP = "Sharp-B.fits" # Similar to flat field, but with sharp edges
BARC = "Arc-B.fits"
RFLAT = "Flat-R.fits" # Median combined flat field for red chip
RSHARP = "Sharp-R.fits" # Similar to flat field, but with sharp edges
RARC = "Arc-R"
XCTEMPLATE = "/Users/cenko/python/lris/std.xcorr.fits"

#######################

def lris_rpreproc(image):

    '''Take an LRIS red image, combine into a single extension, subtract overscan,
    trim, and rotate.'''

    iraf.keck()
    iraf.lris()
    
    # Need grating name to start
    grating = get_head("../rawdata/%s" % image, "GRANAME", extn=0)

    # Convert to single extension fits image
    iraf.multi2simple("../rawdata/%s" % image, "r%s.fits" % image[8:12],
                      overscan=yes, header=yes, trim=yes, verbose=no, debug=no)

    # Trim and rotate
    iraf.imcopy("r%s[%s]" % (image[8:12], LGRATINGS[grating]["trimreg"]),
                "tr%s" % image[8:12])
    iraf.imtranspose("tr%s" % image[8:12], "rtr%s" % image[8:12])

    return

#######################

def lris_bpreproc(image):

    '''Take an LRIS blue image, combine into a single extension, subtract overscan,
    trim, and rotate.'''

    iraf.keck()
    iraf.lris()
    
    # Need grating name to start
    grism = get_head("../rawdata/%s" % image, "GRISNAME", extn=0)

    # Convert to single extension fits image
    iraf.multi2simple("../rawdata/%s" % image, "b%s.fits" % image[8:12],
                      overscan=yes, header=yes, trim=yes, verbose=no, debug=no)

    # Trim and rotate
    iraf.imcopy("b%s[%s]" % (image[8:12], LGRISMS[grism]["trimreg"]),
                "tb%s" % image[8:12])
    iraf.imtranspose("tb%s" % image[8:12], "rtb%s" % image[8:12])

    return

########################

def lris_do_rcalib(flats, arcs):

    '''Calculate y- and x- distortions, create flat-field, calculate wavelength
    solution, for a given grating and central wavelength.'''

    # Pre-process flats
    for image in flats:
        os.system("lris_rpreproc.py %s" % image)

    # Create median combination of flats
    iraf.flatcombine(",".join(["rtr%s.fits" % x[8:12] for x in flats]),
                     output=RFLAT, combine="median", reject="avsigclip",
                     ccdtype="", process=no, subsets=no, delete=no, clobber=no,
                     scale="none", statsec="", lsigma=3.0, hsigma=3.0, rdnoise=0.,
                     gain=1.)

    # Pre-process arcs
    for arc in arcs:
        os.system("lris_rpreproc.py %s" % arc)
        iraf.imcopy("rtr%s.fits" % arc[8:12], "%s%s.fits" % (RARC, arc[8:12]))

    # Generate flats with sharpened edges
    os.system('efits %s "VTKLaplace(i1,numret=1)" %s' % (RFLAT, RSHARP))

    # Calculate y-distortion from sharpened frame
    os.system('getrect -ydist %s -nx 12 -dy 11 -x 1 -y 1 -peakwin 3 -dump' % RSHARP)

    # Copy y rectification to flats
    os.system('copyrect -ydist %s %s' % (RSHARP, RFLAT))
    
    # Find slits
    os.system('findslits -ydist %s -th 0.2 -fwhm 2' % RFLAT)
    
    # Correct Slit locations
    [grating, binning] = get_head(RFLAT, ["GRANAME", "CCDSUM"])
    update_head(RFLAT, "NSLITS", LGRATINGS[grating]["slits"][0])
    update_head(RFLAT, "CSECT1A",
                LGRATINGS[grating]["slits"][1] / int(binning.split()[0]))
    update_head(RFLAT, "CSECT1B",
                LGRATINGS[grating]["slits"][2] / int(binning.split()[0]))
    update_head(RFLAT, "CSECT2A",
                LGRATINGS[grating]["slits"][3] / int(binning.split()[0]))
    update_head(RFLAT, "CSECT2B",
                LGRATINGS[grating]["slits"][4] / int(binning.split()[0]))
                        
    # Create flat fields
    os.system("makeflat -ydist %s" % RFLAT)

    # Correct for bad pixels
    x = pyfits.open("%sflt.fits" % RFLAT.split(".")[0])
    x[0].data[LGRATINGS[grating]["slits"][1]:LGRATINGS[grating]["slits"][2],0:21] = 1
    x[0].data[LGRATINGS[grating]["slits"][3]:LGRATINGS[grating]["slits"][4],0:20] = 1
    x[0].data[LGRATINGS[grating]["slits"][1]:LGRATINGS[grating]["slits"][2],1845:1930] = 1
    x[0].data[LGRATINGS[grating]["slits"][3]:LGRATINGS[grating]["slits"][4],1845:1930] = 1
    x[0].data[LGRATINGS[grating]["slits"][1]:LGRATINGS[grating]["slits"][2],3275:3355] = 1
    x[0].data[LGRATINGS[grating]["slits"][3]:LGRATINGS[grating]["slits"][4],3275:3355] = 1
    x[0].data[LGRATINGS[grating]["slits"][1]:LGRATINGS[grating]["slits"][2],3375:3475] = 1
    x[0].data[LGRATINGS[grating]["slits"][3]:LGRATINGS[grating]["slits"][4],3375:3475] = 1
    x[0].data[LGRATINGS[grating]["slits"][1]:LGRATINGS[grating]["slits"][2],4061:] = 1
    x[0].data[LGRATINGS[grating]["slits"][3]:LGRATINGS[grating]["slits"][4],4061] = 1
    x.writeto("%sflt.fits" % RFLAT.split(".")[0], clobber=True)

    for arc in arcs:

        # Apply distortion, slits and flat-fields to arcs
        os.system("copyrect -ydist %s %s%s.fits" % (RSHARP, RARC, arc[8:12]))
        #os.system("copyslit %s %s -fft -ydist" % (RFLAT, RARC))
        os.system("copyslit %s %s%s.fits -ydist" % (RFLAT, RARC, arc[8:12]))
        #os.system("flat1d %sslt.fits %s -fft -ydist" % (RFLAT.split(".")[0], RARC))
        os.system("flat2d %sflt.fits %s%s.fits" % (RFLAT.split(".")[0], RARC,
                                                   arc[8:12]))

        # Create new files for x-distortion
            #os.system("flat2d %sblz.fits %sf.fits" % (RFLAT.split(".")[0],
        #                                          RARC.split(".")[0]))
        os.system("efits %sblz.fits %s%sf.fits 'i2*greater(i1,0.0)' xdist_R%s.fits" %
                  (RFLAT.split(".")[0], RARC, arc[8:12], arc[8:12]))

        # Calculate x distortion
        os.system("getrect -xdist -ydist xdist_R%s.fits -nx 16 -dy 3 -b 3 -x 2 -y 2 -dump -normwt 1e6" % arc[8:12])

        # Apply x distortion and calculate wavelength solution
        os.system("copyrect -xdist xdist_R%s.fits %s%sf.fits" % (arc[8:12], RARC,
                                                                 arc[8:12]))
        [grating, cwlen] = get_head("../rawdata/%s" % arc, ["GRANAME", "WAVELEN"],
                                    extn=0)
        d0 = LGRATINGS[grating]["dispersion"]
        rlen = get_head("%s%sf.fits" % (RARC, arc[8:12]), "NAXIS1")   
        os.system("waverect -xdist -ydist %s%sf.fits -o 4 -w1 %.2f -w2 %.2f -d0 %.2f -ll %s -arcamp -fwhm 8.5 -bth 20.0" % (RARC, arc[8:12], cwlen - rlen * d0 / 2.0, cwlen + rlen * d0 / 2.0, d0, LRARCLIST))

    return

########################

def lris_do_bcalib(flats, arcs):

    '''Calculate y- and x- distortions, create flat-field, calculate wavelength
    solution, for a given grism.'''

    # Pre-process flats
    for image in flats:
        os.system("lris_bpreproc.py %s" % image)

    # Create median combination of flats
    iraf.flatcombine(",".join(["rtb%s.fits" % x[8:12] for x in flats]),
                     output=BFLAT, combine="median", reject="avsigclip",
                     ccdtype="", process=no, subsets=no, delete=no, clobber=no,
                     scale="none", statsec="", lsigma=3.0, hsigma=3.0, rdnoise=0.,
                     gain=1.)

    # Pre-process arc
    grism = get_head("../rawdata/%s" % arcs[0], "GRISNAME", extn=0)
    for image in arcs:
        os.system("lris_bpreproc.py %s" % image)
    #iraf.imarith("rtb%s.fits" % arcs[0][8:12], "+", "rtb%s.fits" % arcs[1][8:12], BARC)
    os.system("cp rtb%s.fits %s" % (arcs[0][8:12], BARC))
       
    # Generate flats with sharpened edges
    os.system('efits %s "VTKLaplace(i1,numret=1)" %s' % (BFLAT, BSHARP))

    # Calculate y-distortion from sharpened frame
    os.system('getrect -ydist %s -nx 12 -dy 11 -x 1 -y 1 -peakwin 3 -dump' % BSHARP)

    # Copy y rectification to flats
    os.system('copyrect -ydist %s %s' % (BSHARP, BFLAT))
    
    # Find slits
    os.system('findslits -ydist %s -th 0.2 -fwhm 2' % BFLAT)
    
    # Correct Slit locations
    grism = get_head(BFLAT, "GRISNAME")
    update_head(BFLAT, ["NSLITS", "CSECT1A", "CSECT1B", "CSECT2A", "CSECT2B"],
                LGRISMS[grism]["slits"])
                        
    # Create flat fields
    os.system("makeflat -ydist %s" % BFLAT)

    # Bluest part of the flat has too few counts to be useful.
    x = pyfits.open("%sflt.fits" % BFLAT.split(".")[0])
    x[0].data[LGRISMS[grism]["slits"][1]:LGRISMS[grism]["slits"][2],0:1150] = 1
    x[0].data[LGRISMS[grism]["slits"][3]:LGRISMS[grism]["slits"][4],0:1150] = 1
    x.writeto("%sflt.fits" % BFLAT.split(".")[0], clobber=True)

    # Apply distortion, slits and flat-fields to arcs
    os.system("copyrect -ydist %s %s" % (BSHARP, BARC))
    #os.system("copyslit %s %s -fft -ydist" % (BFLAT, BARC))
    os.system("copyslit %s %s -ydist" % (BFLAT, BARC))
    #os.system("flat1d %sslt.fits %s -fft -ydist" % (BFLAT.split(".")[0], BARC))
    os.system("flat2d %sflt.fits %s" % (BFLAT.split(".")[0], BARC))

    # Create new files for x-distortion
    #os.system("flat2d %sblz.fits %sf.fits" % (BFLAT.split(".")[0],
    #                                          BARC.split(".")[0]))
    os.system("efits %sblz.fits %sf.fits 'i2*greater(i1,0.01)' xdist_B.fits" %
              (BFLAT.split(".")[0], BARC.split(".")[0]))

    # Calculate x distortion
    os.system("getrect -xdist -ydist xdist_B.fits -nx 16 -dy 3 -b 3 -x 2 -y 2 -dump -normwt 1e6")

    # Apply x distortion and calculate wavelength solution
    os.system("copyrect -xdist xdist_B.fits %sf.fits" % BARC.split(".")[0])
    d0 = LGRISMS[grism]["dispersion"]
    [w1, w2] = LGRISMS[grism]["wavelengths"]
    os.system("waverect -xdist -ydist %sf.fits -o 5 -w1 %.2f -w2 %.2f -d0 %.2f -ll %s -arcamp -fwhm 8.5 -bth 5.0" % (BARC.split(".")[0], w1, w2, d0, LBARCLIST))

    return

########################

def lris_pipe(rscience, rflats, rarc, bscience, bflats, barcs, docalib=yes,
              doscience=yes):

    '''Take a series of LRIS longslit spectra from a given night, process to
    wavelength-calibrated, sky-subtracted, cosmic-ray cleaned 2D images.'''

    if docalib:
        lris_do_rcalib(rflats, rarc)
        lris_do_bcalib(bflats, barcs)

    if doscience:
        lris_do_rscience(rscience)
        lris_do_bscience(bscience)

    return

##########################

def lris_do_rscience(science):

    '''Using previously created calibration files, process LRIS red longslit
    object (and standard) spectra into clean 2D images.'''
        
    # Pre-process all science images
    scidict = {}
    for image in science:
        [obj, exptime, airmass, skysub, crrej, arc] = get_head("../rawdata/%s" % image,
                                                               ["OBJECT", "EXPTIME",
                                                                "AIRMASS", "SKYSUB",
                                                                "CRREJ", "ARCNAME"],
                                                               extn=0)
        if not (skysub == 1):
            skysub = 0
        if not (crrej == 1):
            crrej = 0
        obj = obj.replace(" ", "_")
        if not scidict.has_key(obj):
            scidict[obj] = [[image, exptime, airmass, skysub, crrej, arc]]
        else:
            scidict[obj].append([image, exptime, airmass, skysub, crrej, arc])
        os.system("lris_rpreproc.py %s" % image)

    # Now process all science frame
    for obj, dets in scidict.iteritems():

        for i in range(len(dets)):

            [image, exptime, airmass, skysub, crrej, narc] = dets[i]
            
            rsci = "rtr%s.fits" % image[8:12]
        
            # Apply y distortion
            os.system("copyrect -ydist %s %s" % (RSHARP, rsci))
        
            # Copy slits
            os.system("copyslit %s %s -ydist" % (RFLAT, rsci))
        
            # Apply flats
            #os.system("flat1d %sslt.fits %s -fft -ydist" % (RFLAT.split(".")[0], rsci))
            os.system("flat2d %sflt.fits %s" % (RFLAT.split(".")[0], rsci))
            #os.system("flat2d %sblz.fits %sf.fits" % (RFLAT.split(".")[0],
                                                      #rsci.split(".")[0]))

            # Apply x distortion
            os.system("copyrect -xdist xdist_R%s.fits %sf.fits" % (narc,
                                                                   rsci.split(".")[0]))

            # Apply wavelength solution
            os.system("copyrect -wdist %s%sf.fits %sf.fits" % (RARC, narc,
                                                               rsci.split(".")[0]))

            # Tweak based on sky lines (if possible)
            os.system("waverect -wdist -ydist %sf.fits -ll %s -zero 1 -skyamp -bth 5.0 -fwhm 8.5 -slit 2" % (rsci.split(".")[0], LSKYLIST))

            # If desired, subtract sky
            os.system("skyrect -ydist %sf.fits -skyint 512 -kx 3 -ky 5 -b 5 -dkx 0.75 -skycut 5.0 -check 15 -pct 40 -g %.2f -rn %.2f" % (rsci.split(".")[0], RCCDGAIN, RCCDRNOISE))
            if skysub:
                os.system("efits %sf.fits %sfm.fits 'i1-i2' %sfo.fits" %
                          (rsci.split(".")[0], rsci.split(".")[0], rsci.split(".")[0]))
            else:
                shutil.copy("%sf.fits" % rsci.split(".")[0],
                            "%sfo.fits" % rsci.split(".")[0])

            # Clean cosmic rays
            if crrej:
                os.system("pyrclean %sfo.fits -th 7.0 -rn %f -g %f -M 65000.0 -bigsky" % (rsci.split(".")[0], RCCDRNOISE, RCCDGAIN))
            else:
                shutil.copy("%sfo.fits" % rsci.split(".")[0],
                            "%sfoc.fits" % rsci.split(".")[0])

            # Unrectify
            os.system("unrect -wdist -xdist -ydist %sfoc.fits" % rsci.split(".")[0])

            # Rename and add / remove keywords
            shutil.copy("%sfocr.fits" % rsci.split(".")[0],
                        "%s_%02i_R.fits" % (obj, i+1))
            update_head("%s_%02i_R.fits" % (obj, i+1), ["EXPTIME", "AIRMASS", "CRREJ",
                                                        "SKYSUB", "ONAME", "DISPAXIS",
                                                        "OBJECT", "OBSERVAT"],
                                                       [exptime, airmass, crrej,
                                                        skysub, image, 1, obj, "Keck"])
            delete_head("%s_%02i_R.fits" % (obj, i+1), ["LTM1_2", "LTM2_1",
                                                           "CD1_2", "CD2_1"])

    return

#################################################

def lris_do_bscience(science):

    '''Using previously created calibration files, process LRIS blue longslit
    object (and standard) spectra into clean 2D images.'''
        
    # Pre-process all science images
    scidict = {}
    for image in science:
        update_head("../rawdata/%s" % image, "EXPTIME",
                    get_head("../rawdata/%s" % image, "TTIME", extn=0))
        [obj, exptime, airmass, skysub, crrej] = get_head("../rawdata/%s" % image,
                                                          ["OBJECT", "EXPTIME",
                                                           "AIRMASS", "SKYSUB",
                                                           "CRREJ"], extn=0)
        if not (skysub == 1):
            skysub = 0
        if not (crrej == 1):
            crrej = 0
        obj = obj.replace(" ", "_")
        if not scidict.has_key(obj):
            scidict[obj] = [[image, exptime, airmass, skysub, crrej]]
        else:
            scidict[obj].append([image, exptime, airmass, skysub, crrej])
        os.system("lris_bpreproc.py %s" % image)

    # Now process all science frame
    for obj, dets in scidict.iteritems():

        for i in range(len(dets)):

            [image, exptime, airmass, skysub, crrej] = dets[i]
            
            bsci = "rtb%s.fits" % image[8:12]
        
            # Apply y distortion
            os.system("copyrect -ydist %s %s" % (BSHARP, bsci))
        
            # Copy slits
            os.system("copyslit %s %s -ydist" % (BFLAT, bsci))
        
            # Apply flats
            #os.system("flat1d %sslt.fits %s -fft -ydist" % (BFLAT.split(".")[0], bsci))
            os.system("flat2d %sflt.fits %s" % (BFLAT.split(".")[0], bsci))
            #os.system("flat2d %sblz.fits %sf.fits" % (BFLAT.split(".")[0],
                                                      #bsci.split(".")[0]))

            # Apply x distortion
            os.system("copyrect -xdist xdist_B.fits %sf.fits" % bsci.split(".")[0])

            # Apply wavelength solution
            os.system("copyrect -wdist %sf.fits %sf.fits" % (BARC.split(".")[0],
                                                              bsci.split(".")[0]))

            # Tweak based on sky lines (if possible)
            #os.system("waverect -wdist -ydist %sf.fits -ll %s -zero 1 -skyamp -bth 1.5" % (bsci.split(".")[0], LSKYLIST))

            # If desired, subtract sky
            if skysub:
                os.system("skyrect -ydist %sf.fits -skyint 512 -kx 3 -ky 1 -b 5 -dkx 0.75 -skycut 5.0 -check 5 -pct 40" % bsci.split(".")[0])
                os.system("efits %sf.fits %sfm.fits 'i1-i2' %sfo.fits" %
                          (bsci.split(".")[0], bsci.split(".")[0], bsci.split(".")[0]))
            else:
                shutil.copy("%sf.fits" % bsci.split(".")[0],
                            "%sfo.fits" % bsci.split(".")[0])

            # Clean cosmic rays
            if crrej:
                os.system("pyrclean %sfo.fits -th 7.0 -rn %f -g %f -M 65000.0 -bigsky -k 1 -mx 15 -my 15" % (bsci.split(".")[0], BCCDRNOISE, BCCDGAIN))
            else:
                shutil.copy("%sfo.fits" % bsci.split(".")[0],
                            "%sfoc.fits" % bsci.split(".")[0])

            # Unrectify
            os.system("unrect -wdist -xdist -ydist %sfoc.fits" % bsci.split(".")[0])

            # Rename and add / remove keywords
            shutil.copy("%sfocr.fits" % bsci.split(".")[0],
                        "%s_%02i_B.fits" % (obj, i+1))
            update_head("%s_%02i_B.fits" % (obj, i+1), ["EXPTIME", "AIRMASS", "CRREJ",
                                                        "SKYSUB", "ONAME", "DISPAXIS",
                                                        "OBJECT", "OBSERVAT"],
                                                       [exptime, airmass, crrej,
                                                        skysub, image, 1, obj, "Keck"])
            delete_head("%s_%02i_B.fits" % (obj, i+1), ["LTM1_2", "LTM2_1",
                                                           "CD1_2", "CD2_1"])

    return

#################################################

def lris_extract(science, standards, dostandard=yes):

    '''Extract and flux calibrate LRIS spectra (assuming they have been
    processed into 2D images using LRIS_pipe above).'''

    if dostandard:
        lris_standard(standards)

    for src in science:

        bimages = iraffiles("%s_??_B.fits" % src)
        rimages = iraffiles("%s_??_R.fits" % src)
        joinstr = ""

        for bimage in bimages:

            bimage = bimage.split(".")[0]
            
            # Find appropriate standard for this image
            bstd = get_head("%s.fits" % bimage, "STDNAME")
            
            # Extract 1D spectra
            if get_head("%s.fits" % bimage, "SKYSUB"):
                iraf.apall(bimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=no, background="none",
                           reference=bstd)
            else:
                iraf.apall(bimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=no, background="fit",
                           reference=bstd)

            # Normalize by the standard continuua
            iraf.sarith("%s.ms" % bimage, "/", "../standards/s%s.ms.fits" % bstd,
                        "s%s.ms" % bimage, w1=INDEF, w2=INDEF)            
            
            # Telluric correct
            iraf.telluric("s%s.ms" % bimage, "ts%s.ms" % bimage,
                          "../standards/telluric.%s.fits" % bstd, tweakrms=yes,
                          interactive=yes, sample='6850:6950,7575:7700,8750:9925')
            
            # Flux calibration
            iraf.calibrate("ts%s.ms" % bimage, "fts%s.ms" % bimage, extinct=yes,
                           flux=yes, extinction="home$extinct/maunakeaextinct.dat",
                           observatory="Keck", sens="../standards/%s.sens" % bstd,
                           airmass='', exptime='')

            joinstr += "fts%s.ms," % bimage

        for rimage in rimages:

            rimage = rimage.split(".")[0]
            
            # Find appropriate standard for this image
            rstd = get_head("%s.fits" % rimage, "STDNAME")
            
            # Extract 1D spectra
            if get_head("%s.fits" % rimage, "SKYSUB"):
                iraf.apall(rimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=no, background="none",
                           reference=rstd)
            else:
                iraf.apall(rimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=no, background="fit",
                           reference=rstd)

            # Normalize by the standard continuua
            iraf.sarith("%s.ms" % rimage, "/", "../standards/s%s.ms.fits" % rstd,
                        "s%s.ms" % rimage, w1=5500, w2=INDEF)
            
            # Telluric correct
            iraf.telluric("s%s.ms" % rimage, "ts%s.ms" % rimage,
                          "../standards/telluric.%s.fits" % rstd, tweakrms=yes,
                          interactive=yes, sample='6850:6950,7575:7700')
            
            # Flux calibration
            iraf.calibrate("ts%s.ms" % rimage, "fts%s.ms" % rimage, extinct=yes,
                           flux=yes, extinction="home$extinct/maunakeaextinct.dat",
                           observatory="Keck", sens="../standards/%s.sens" % rstd,
                           airmass='', exptime='')

            joinstr += "fts%s.ms," % rimage
           
        # Combine
        iraf.scombine(joinstr[:-1], "%s.ms.fits" % src, combine="average",
                      reject="avsigclip", w1=INDEF, w2=INDEF, dw=INDEF, nw=INDEF,
                      scale="median", zero="none", weight="none", sample="5450:5600",
                      lsigma=3.0, hsigma=3.0, gain=RCCDGAIN, rdnoise=RCCDRNOISE)

        # Plot to enable final tweaks
        iraf.splot("%s.ms.fits" % src)

    return

#################################################

def lris_standard(standards, xcorr=yes):

    '''Extract standard stars and calculate sensitivity functions.'''

    for ofile, nstd in standards:

        shutil.copy(ofile, "%s.fits" % nstd)
        
        # Extract standard
        if get_head("%s.fits" % nstd, "SKYSUB"):
            iraf.apall(nstd, output="", inter=yes, find=yes, recenter=yes, resize=yes,
                       edit=yes, trace=yes, fittrace=yes, extract=yes, extras=yes,
                       review=no, background="none")
        else:
            iraf.apall(nstd, output="", inter=yes, find=yes, recenter=yes, resize=yes,
                       edit=yes, trace=yes, fittrace=yes, extract=yes, extras=yes,
                       review=no, background="fit")

        if xcorr:
            # Cross correlate to tweak wavelength
            iraf.xcsao("%s.ms" % nstd, templates=XCTEMPLATE,
                       correlate="wavelength", logfiles="%s.xcsao" % ofile)

            # If a solution was found, apply shift
            try:
                xcsaof = open("%s.xcsao" % ofile)
                shift = float(xcsaof.readlines()[6].split()[14])
            except:
                shift = 0.0
                iraf.specshift("%s.ms" % nstd, -shift)
        
        # Create telluric
        iraf.splot("%s.ms" % nstd) # Remove telluric and absorption, save as a%s.ms
        iraf.sarith("%s.ms" % nstd, "/", "a%s.ms" % nstd, "telluric.%s.fits" % nstd)
        iraf.imreplace("telluric.%s.fits" % nstd, 1.0, lower=0.0, upper=0.0)
        iraf.splot("telluric.%s.fits" % nstd) # Remove stellar features and resave
        
        # Create smoothed standard
        iraf.gauss("a%s.ms[*,1,1]" % nstd, "s%s.ms" % nstd, 5.0)
        iraf.sarith("%s.ms" % nstd, "/", "s%s.ms" % nstd, "ds%s.ms" % nstd)

        # Apply telluric correction
        iraf.telluric("ds%s.ms" % nstd, "tds%s.ms" % nstd, "telluric.%s.fits" % nstd,
                      xcorr=no, tweakrms=no, interactive=no,
                      sample='4000:4010,6850:6975,7150:7350,7575:7725,8050:8400,8900:9725')

        # Define bandpasses for standard star calculation
        obj=get_head("%s.fits" % nstd, "OBJECT")
        iraf.standard("tds%s.ms" % nstd, "%s.std" % nstd,
                  extinction='home$extinct/maunakeaextinct.dat',
                  caldir='home$standards/', observatory='Keck', interac=yes,
                  star_name=obj, airmass='', exptime='')
    
        # Determine sensitivity function
        iraf.sensfunc('%s.std' % nstd, '%s.sens' % nstd,
                      extinction='home$extinct/maunakeaextinct.dat',
                      newextinction='extinct.dat', observatory='Keck',
                      function='legendre', order=4, interactive=yes)

    return

##########################

def lris_rename(oldname, newname):

    '''Rename all the files associated with a given target'''

    inlis = iraffiles("*%s*" % oldname)

    for f in inlis:
        newf = f.replace(oldname, newname)
        os.system("mv %s %s" % (f, newf))

############################

def lris_plot(flfile):

    import numpy as np
    import matplotlib.pylab as plt

    x = np.loadtxt(flfile)
    plt.clf()
    plt.plot(x[:,0], x[:,1])
    plt.xlabel("Observed Wavelength ($\AA$)", fontsize=20)
    plt.ylabel("$f_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", fontsize=20)
    plt.title(flfile.split("-")[0], fontsize=24)
    plt.xlim(np.min(x[:,0])-100.0, np.max(x[:,0]+100.0))
    plt.ylim(0.0, 1.1*np.max(x[:,1]))
    plt.savefig("%s.png" % flfile[:-4])
    plt.show()










